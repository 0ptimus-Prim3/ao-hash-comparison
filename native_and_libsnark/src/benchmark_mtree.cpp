#include "gadget/hash/anemoi/anemoi_gadget.hpp"
#include "gadget/hash/arion/arion_gadget.hpp"
#include "gadget/hash/griffin/griffin_gadget.hpp"
#include "gadget/hash/mimc/mimc256_gadget.hpp"
#include "gadget/hash/mimc/mimc512f2k_gadget.hpp"
#include "gadget/hash/mimc/mimc512f_gadget.hpp"
#include "gadget/hash/poseidon/poseidon_gadget.hpp"
#include "gadget/hash/poseidon2/poseidon2_gadget.hpp"
#include "gadget/hash/rescue/rescue_gadget.hpp"
#include "gadget/hash/sha256/sha256_gadget_pp.hpp"
#include "gadget/hash/sha512/sha512_gadget_pp.hpp"
#include "gadget/mtree_gadget.hpp"

#include "hash/anemoi/anemoi.hpp"
#include "hash/arion/arion.hpp"
#include "hash/griffin/griffin.hpp"
#include "hash/poseidon/poseidon.hpp"
#include "hash/poseidon2/poseidon2.hpp"
#include "hash/rescue/rescue.hpp"

#include "r1cs/r1cs_ppzksnark_pp.hpp"
#include "tree/mtree.hpp"
#include "util/measure.hpp"

#include <filesystem>
#include <fstream>
#include <libff/common/default_types/ec_pp.hpp>
#include <libsnark/common/default_types/r1cs_ppzksnark_pp.hpp>
#include <libsnark/zk_proof_systems/ppzksnark/r1cs_ppzksnark/r1cs_ppzksnark.hpp>
#include <string>

static constexpr size_t MIN_HEIGHT = 6;
static constexpr size_t MAX_HEIGHT = 30 + 1; // the +1 is to highlight that the bound is exclusive
static constexpr size_t STEP_HEIGHT = 6;
static constexpr int NUM_THREADS = 1;

namespace fs = std::filesystem;

using ppT = libsnark::default_r1cs_ppzksnark_pp;
using FieldT = libff::Fr<ppT>;

static std::ofstream log_file;
static std::string table_header = std::string("Height\t") + //std::string("Tree\t") +
                                  //std::string("Gadget\t") + std::string("Constraint\t") +
                                  //std::string("Witness\t") + std::string("Key\t") +
                                  std::string("Proof\n"); // + std::string("Verify\n");

template<size_t height, typename GadHash>
bool bench_mtree(size_t trans_idx = 0)
{
    static constexpr size_t HEIGHT = height;

    using GadTree = MTreeGadget<height, GadHash>;
    using DigVar = typename GadTree::DigVar;
    using Level = typename GadTree::Level;
    using Hash = typename GadHash::Hash;
    using Tree = MTreePath<HEIGHT, Hash>;

    static constexpr size_t DIGEST_VARS = GadHash::DIGEST_VARS;

    static std::mt19937 rng{std::random_device{}()};

    double elap = 0;
    std::vector<uint8_t> data(Tree::INPUT_SIZE);
    std::generate(data.begin(), data.end(), std::ref(rng));
    Tree tree;

    log_file << height << '\t';
    log_file.flush();

    // Build tree
    elap = measure(
        [&]() {
            tree = Tree{data.begin(), data.end()};
        },
        1, 1, "Tree Generation", false);
    //log_file << elap << '\t';
    log_file.flush();

    // Test Gadget
    libsnark::protoboard<FieldT> pb;

    DigVar out{pb, DIGEST_VARS, FMT("out")};
    DigVar trans{pb, DIGEST_VARS, FMT("trans")};
    std::vector<Level> other;
    std::unique_ptr<GadTree> gadget;

    for (size_t i = 0; i < HEIGHT - 1; ++i)
        other.emplace_back(make_uniform_array<Level>(pb, DIGEST_VARS, FMT("other_%llu", i)));

    // Gadget construction
    elap = measure(
        [&]() { gadget = std::make_unique<GadTree>(pb, out, trans, other, FMT("merkle_tree")); }, 1,
        1, "Gadget construction", false);
    //log_file << elap << '\t';
    log_file.flush();

    pb.set_input_sizes(DIGEST_VARS);

    // Constraint generation
    elap = measure(
        [&]()
        {
            out.generate_r1cs_constraints();
            trans.generate_r1cs_constraints();
            for (size_t i = 0; i < other.size(); ++i)
                for (size_t j = 0; j < other[i].size(); ++j)
                    other[i][j].generate_r1cs_constraints();
            gadget->generate_r1cs_constraints();
        },
        1, 1, "Constraint generation", false);
    //log_file << elap << '\t';
    log_file.flush();

    trans.generate_r1cs_witness(tree.get_node(0)->get_digest());

    for (size_t i = 0, idx = trans_idx; i < other.size(); ++i, idx /= Tree::ARITY)
    {
        size_t j = idx % Tree::ARITY;
        size_t off = Hash::BLOCK_SIZE + i * (Tree::ARITY - 1) * Hash::DIGEST_SIZE;

        for (size_t k = 0; k < j; ++k)
            other[i][k].generate_r1cs_witness(data.data() + off + k * Hash::DIGEST_SIZE,
                                              Hash::DIGEST_SIZE);

        other[i][j].generate_r1cs_witness(tree.get_node(i)->get_digest());

        for (size_t k = j + 1; k < Tree::ARITY; ++k)
            other[i][k].generate_r1cs_witness(data.data() + off + (k - 1) * Hash::DIGEST_SIZE,
                                              Hash::DIGEST_SIZE);
    }

    // Witness generation
    elap = measure([&]() { gadget->generate_r1cs_witness(trans_idx); }, 1, 1, "Witness generation",
                   false);
    //log_file << elap << '\t';
    log_file.flush();

    // Key generation
    r1cs_ppzksnark_keypair<ppT> keypair;
    elap = measure(
        [&]() { keypair = libsnark::r1cs_ppzksnark_generator<ppT>(pb.get_constraint_system()); }, 1,
        1, "Key generation", false);
    //log_file << elap << '\t';
    log_file.flush();

    // Proof generation
    libsnark::r1cs_ppzksnark_proof<ppT> proof;
    elap = measure(
        [&]()
        {
            proof = libsnark::r1cs_ppzksnark_prover<ppT>(keypair.pk, pb.primary_input(),
                                                         pb.auxiliary_input());
        },
        1, 4, "Proof generation", false);
    log_file << elap << '\n';
    log_file.flush();

    // Proof Verification
    bool result;
    elap = measure(
        [&]()
        {
            result = libsnark::r1cs_ppzksnark_verifier_strong_IC<ppT>(keypair.vk,
                                                                      pb.primary_input(), proof);
        },
        1, 1, "Proof verification", false);
    //log_file << elap << '\n';
    log_file.flush();


    return result;
}

constexpr size_t countr_zero(size_t x)
{
    size_t result = 0;
    while ((x & 1) == 0)
    {
        x >>= 1;
        ++result;
    }

    return result;
}

template<size_t first, size_t last, typename GadHash>
void bench_range(const char *name)
{
    static constexpr size_t ARITY = GadHash::BLOCK_SIZE / GadHash::DIGEST_SIZE;

    if constexpr (first < last)
    {
        bench_mtree<first / countr_zero(ARITY), GadHash>();

        bench_range<first + STEP_HEIGHT, last, GadHash>(name);
    }
}

template<typename GadHash>
void bench(const char *name)
{
    static constexpr size_t RATIO = GadHash::BLOCK_SIZE / GadHash::DIGEST_SIZE;

    log_file << name << " (" << RATIO << ":1), r = " << GadHash::Hash::ROUNDS_N
             << ", c = " << GadHash::size() << '\n';
    log_file << table_header;
    bench_range<MIN_HEIGHT, MAX_HEIGHT, GadHash>(name);
    log_file << '\n';
}

int main()
{
    fs::create_directories("./log");
    std::string timestamp = std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(
                                               std::chrono::system_clock::now().time_since_epoch())
                                               .count());
    std::string log_file_name = std::string("./log/benchmark_mtree_") + timestamp +
                                std::string(".log");
    log_file.open(log_file_name);

    std::cout << "Logging to " << log_file_name << "...\n";

    log_file << std::boolalpha;
    libff::inhibit_profiling_info = true;
    libff::inhibit_profiling_counters = true;

    libff::default_ec_pp::init_public_params();

    log_file << "Merkle Tree Benchmark"
             << "\n";
    log_file << "Prime:\t" << FieldT::mod << "\n";
    log_file << "d:\t"
             << "5"
             << "\n";
    log_file << "Minimum Merkle Tree Height:\t" << MIN_HEIGHT << "\n";
    log_file << "Maximum Merkle Tree Height:\t" << MAX_HEIGHT - 1 << "\n\n";

#ifdef MULTICORE
    omp_set_num_threads(NUM_THREADS);
    log_file << "Threads:\t" << omp_get_max_threads() << "\n";
#else
    log_file << "Threads:\t1\n";
#endif

    /**
    To benchmark some permutation gadget over a merkle tree, follow the syntax below:

    log_file << "Permutation_Name\n";
    log_file << table_header;
    bench_range<min_tree_height, max_tree_height, PermutationGadget>("Permutation_Name"); 
    lof_file << "\n";

    Different gadgets can require different template arguments to be instantiated. Refer to their 
    documentation for more details. 
    Example:
    GriffinGadget, PoseidonGadget and ArionGadget require an instantiation of the respective
    permutation as argument. 
    In particular, an instantiation of ArionGadget follows the syntax below:
        PermutationGadget = ArionGadget<Arion<FieldT, rate, capacity, rounds>>
    Griffin follows the same syntax, while Poseidon require both the number of full rounds and the 
    number of partial rounds.
    **/

    bench<ArionGadget<Arion<FieldT, 2, 1, 6>>>("Arion");
    bench<ArionGadget<Arion<FieldT, 4, 1, 5>>>("Arion");
    bench<ArionGadget<Arion<FieldT, 8, 1, 4>>>("Arion");

    bench<AnemoiGadget<Anemoi<FieldT, 2, 2, 14>>>("Anemoi");
    bench<AnemoiGadget<Anemoi<FieldT, 4, 2, 12>>>("Anemoi");
    bench<AnemoiGadget<Anemoi<FieldT, 8, 2, 11>>>("Anemoi");

    bench<GriffinGadget<Griffin<FieldT, 2, 1, 12>>>("Griffin");
    bench<GriffinGadget<Griffin<FieldT, 4, 4, 9>>>("Griffin");
    bench<GriffinGadget<Griffin<FieldT, 8, 4, 9>>>("Griffin");

    bench<PoseidonGadget<Poseidon<FieldT, 2, 1, 4, 55>>>("Poseidon");
    bench<PoseidonGadget<Poseidon<FieldT, 4, 1, 4, 56>>>("Poseidon");
    bench<PoseidonGadget<Poseidon<FieldT, 8, 1, 4, 56>>>("Poseidon");

    bench<Poseidon2Gadget<Poseidon2<FieldT, 2, 4, 55>>>("Poseidon2");
    bench<Poseidon2Gadget<Poseidon2<FieldT, 4, 4, 56>>>("Poseidon2");
    bench<Poseidon2Gadget<Poseidon2<FieldT, 8, 4, 56>>>("Poseidon2");

    bench<RescueGadget<Rescue<FieldT, 2, 1, 14>>>("Rescue");
    bench<RescueGadget<Rescue<FieldT, 4, 1, 9>>>("Rescue");
    bench<RescueGadget<Rescue<FieldT, 8, 1, 8>>>("Rescue");

    bench<Sha256Gadget<FieldT>>("SHA-256");

    log_file.close();

    return 0;
}
