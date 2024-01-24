#include "gadget/mtree_gadget.hpp"
#include "gadget/hash/arion/arion_gadget.hpp"
#include "gadget/hash/griffin/griffin_gadget.hpp"
#include "gadget/hash/mimc/mimc256_gadget.hpp"
#include "gadget/hash/mimc/mimc512f2k_gadget.hpp"
#include "gadget/hash/mimc/mimc512f_gadget.hpp"
#include "gadget/hash/poseidon/poseidon_gadget.hpp"
#include "gadget/hash/sha256/sha256_gadget_pp.hpp"
#include "gadget/hash/sha512/sha512_gadget_pp.hpp"

#include "hash/anemoi/anemoi.hpp"
#include "hash/arion/arion.hpp"
#include "hash/griffin/griffin.hpp"
#include "hash/poseidon/poseidon.hpp"

#include "tree/mtree.hpp"
#include "util/measure.hpp"

#include <libsnark/common/default_types/r1cs_ppzksnark_pp.hpp>
#include <libsnark/zk_proof_systems/ppzksnark/r1cs_ppzksnark/r1cs_ppzksnark.hpp>

#include <libff/common/default_types/ec_pp.hpp>

#include <fstream>
#include <omp.h>

using ppT = libsnark::default_r1cs_ppzksnark_pp;
using FieldT = libff::Fr<ppT>;

template<typename GadTree, bool full_tree = false>
bool test_mtree()
{
    static constexpr size_t HEIGHT = GadTree::HEIGHT;

    using DigVar = typename GadTree::DigVar;
    using Level = typename GadTree::Level;
    using GadHash = typename GadTree::GadHash;
    using Hash = typename GadHash::Hash;
    using Tree = std::conditional_t<full_tree, MTree<HEIGHT, Hash>, MTreePath<HEIGHT, Hash>>;
    using Node = typename Tree::Node;

    static constexpr size_t DIGEST_VARS = GadHash::DIGEST_VARS;

    static std::mt19937 rng{std::random_device{}()};

    // Build tree
    std::vector<uint8_t> data(Tree::INPUT_SIZE);
    std::generate(data.begin(), data.end(), std::ref(rng));

    if constexpr (!GadTree::HASH_ISBOOLEAN)
        field_clamp<FieldT>(data.data(), data.size());

    size_t trans_idx = 0;
    Tree tree{[&]()
              {
                  if constexpr (full_tree)
                      return Tree{data.begin(), data.end()};
                  else
                      return Tree{data.begin(), data.end(), trans_idx};
              }()};

    std::cout << '\n' << tree;

    // There is no need to extract as FieldVariable does that for us

    // Test Gadget
    libsnark::protoboard<FieldT> pb;

    DigVar out{pb, DIGEST_VARS, FMT("out")};
    DigVar trans{pb, DIGEST_VARS, FMT("trans")};
    std::vector<Level> other;

    for (size_t i = 0; i < HEIGHT - 1; ++i)
        other.emplace_back(make_uniform_array<Level>(pb, DIGEST_VARS, FMT("other_%llu", i)));

    GadTree gadget{pb, out, trans, other, FMT("merkle_tree")};

    pb.set_input_sizes(DIGEST_VARS);
    out.generate_r1cs_constraints();
    trans.generate_r1cs_constraints();
    for (size_t i = 0; i < other.size(); ++i)
        for (size_t j = 0; j < other[i].size(); ++j)
            other[i][j].generate_r1cs_constraints();
    gadget.generate_r1cs_constraints();

    if constexpr (full_tree)
    {
        trans.generate_r1cs_witness(tree.get_node(trans_idx)->get_digest());

        const Node *aux = tree.get_node(trans_idx)->parent();
        for (size_t i = 0; i < other.size(); ++i, aux = aux->parent())
            for (size_t j = 0; j < other[i].size(); ++j)
                other[i][j].generate_r1cs_witness(aux->child(j)->get_digest());
    }
    else
    {
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
    }

    gadget.generate_r1cs_witness(trans_idx);

    std::string vanilla_dump{hexdump(tree.digest(), Hash::DIGEST_SIZE)};
    std::string zkp_dump;

    if constexpr (GadTree::HASH_ISBOOLEAN)
    {
        uint8_t buff[Hash::DIGEST_SIZE]{};
        std::vector<bool> buff_bv(DIGEST_VARS);

        for (size_t i = 0; i < DIGEST_VARS; ++i)
            buff_bv[i] = pb.val(out.bits[i]).as_ulong();
        pack_bits(buff, buff_bv);

        zkp_dump = hexdump(buff);
    }
    else
    {

        for (auto &&x : out)
        {
            zkp_dump += hexdump(pb.val(x));
        }
    }

    std::cout << "\nTrans idx: " << trans_idx << '\n';
    std::cout << "Vanilla output:\t" << vanilla_dump << '\n';
    std::cout << "ZKP output:\t" << zkp_dump << '\n';

    bool result = vanilla_dump == zkp_dump;

    auto keypair{libsnark::r1cs_ppzksnark_generator<ppT>(pb.get_constraint_system())};
    auto proof{
        libsnark::r1cs_ppzksnark_prover<ppT>(keypair.pk, pb.primary_input(), pb.auxiliary_input())};

    result &= libsnark::r1cs_ppzksnark_verifier_strong_IC<ppT>(keypair.vk, pb.primary_input(),
                                                               proof);

    return result;
}

static bool run_tests()
{
    static constexpr size_t TREE_HEIGHT = 4;

    bool check = true;
    bool all_check = true;
    std::cout << std::boolalpha;
    libff::inhibit_profiling_info = true;
    libff::inhibit_profiling_counters = true;

    ppT::init_public_params();
    /*
    std::cout << "SHA256... ";
    std::cout.flush();
    {
        check = test_mtree<MTreeGadget<TREE_HEIGHT, Sha256Gadget<FieldT>>, true>();
    }
    std::cout << check << '\n';
    all_check &= check;

    std::cout << "Path SHA256... ";
    std::cout.flush();
    {
        check = test_mtree<MTreeGadget<TREE_HEIGHT, Sha256Gadget<FieldT>>>();
    }
    std::cout << check << '\n';
    all_check &= check;

    std::cout << "MiMC256... ";
    std::cout.flush();
    {
        //      check = test_pmtree<TREE_HEIGHT, GadMimc256>();
    }
    std::cout << check << '\n';
    all_check &= check;

    std::cout << "Path MiMC256... ";
    std::cout.flush();
    {
        //     check = test_pmtree_path<TREE_HEIGHT, GadMimc256>();
    }
    std::cout << check << '\n';
    all_check &= check;

    std::cout << "Path SHA512... ";
    std::cout.flush();
    {
        check = test_mtree<MTreeGadget<TREE_HEIGHT, Sha512Gadget<FieldT>>>();
    }
    std::cout << check << '\n';
    all_check &= check;

    std::cout << "Path MiMC512F... ";
    std::cout.flush();
    {
        //      check = test_pmtree_path<TREE_HEIGHT, GadMimc512F>();
    }
    std::cout << check << '\n';
    all_check &= check;

    std::cout << "Path MiMC512F2K... ";
    std::cout.flush();
    {
        //      check = test_pmtree_path<TREE_HEIGHT, GadMimc512F2K>();
    }
    std::cout << check << '\n';
    all_check &= check;

    std::cout << "Path Griffin... ";
    std::cout.flush();
    {
        check = test_mtree<MTreeGadget<TREE_HEIGHT, GriffinGadget<Griffin<FieldT, 2, 1>>>>();
    }
    std::cout << check << '\n';
    all_check &= check;


    std::cout << "Path Poseidon... ";
    std::cout.flush();
    {
        check = test_mtree<MTreeGadget<TREE_HEIGHT, PoseidonGadget<Poseidon<FieldT, 2, 1>>>>();
    }
    std::cout << check << '\n';
    all_check &= check;*/

    std::cout << "Arion... ";
    std::cout.flush();
    {
        check = test_mtree<MTreeGadget<TREE_HEIGHT, ArionGadget<Arion<FieldT, 4, 1>>>, true>();
    }
    std::cout << check << '\n';
    all_check &= check;

    std::cout << "Path Arion... ";
    std::cout.flush();
    {
        check = test_mtree<MTreeGadget<TREE_HEIGHT, ArionGadget<Arion<FieldT, 4, 1>>>, false>();
    }
    std::cout << check << '\n';
    all_check &= check;


    return all_check;
}

int main()
{
    std::cout << "\n==== Testing MerkleTree Gadget ====\n";

    bool all_check = run_tests();

    std::cout << "\n==== " << (all_check ? "ALL TESTS SUCCEEDED" : "SOME TESTS FAILED")
              << " ====\n\n";

#ifdef MEASURE_PERFORMANCE
#endif

    return 0;
}
