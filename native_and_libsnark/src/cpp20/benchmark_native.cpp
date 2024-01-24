#include "hash/anemoi/anemoi.hpp"
#include "hash/arion/arion.hpp"
#include "hash/griffin/griffin.hpp"
#include "hash/poseidon/poseidon.hpp"
#include "hash/poseidon2/poseidon2.hpp"
#include "hash/rescue/rescue.hpp"
#include "hash/sha/sha256.hpp"
#include "tree/mtree.hpp"
#include "util/measure.hpp"

#include <chrono>
#include <filesystem>
#include <fstream>
#include <omp.h>
#include <string>


static constexpr size_t MIN_HEIGHT = 6;
static constexpr size_t MAX_HEIGHT = 18 + 1; // the +1 is to highlight that the bound is exclusive
static constexpr size_t STEP_HEIGHT = 6;
static constexpr int NUM_THREADS = 1;

namespace fs = std::filesystem;

using FieldT = libff::Fr<libff::default_ec_pp>;

static std::ofstream log_file;
static std::string table_header = std::string("Height\t") + std::string("Time");

template<size_t height, typename Hash>
void bench_mtree(size_t trans_idx = 0)
{
    static constexpr size_t HEIGHT = height;

    using Tree = MTree<HEIGHT, Hash>;

    std::mt19937 rng{std::random_device{}()};
    double elap = 0;
    std::unique_ptr<Tree> tree;
    std::vector<uint8_t> data(Tree::INPUT_SIZE);
    //std::generate(data.begin(), data.end(), std::ref(rng));

    log_file << height << '\t';
    log_file.flush();

    // Build tree
    elap = measure([&]() { tree = std::make_unique<Tree>(data.begin(), data.end()); }, 1, 1,
                   "Tree Generation", false);

    log_file << elap << '\n';
    log_file.flush();
}

template<size_t first, size_t last, typename Hash>
void bench_range(const char *name)
{
    static constexpr size_t ARITY = Hash::BLOCK_SIZE / Hash::DIGEST_SIZE;

    if constexpr (first < last)
    {
        bench_mtree<first / std::countr_zero(ARITY), Hash>();

        bench_range<first + STEP_HEIGHT, last, Hash>(name);
    }
}

template<typename Hash>
void bench(const char *name)
{
    static constexpr size_t RATIO = Hash::BLOCK_SIZE / Hash::DIGEST_SIZE;
    static const char *KIND_NAME = "";

    log_file << name << ' ' << KIND_NAME << " (" << RATIO << ":1), r = " << Hash::ROUNDS_N << '\n';
    log_file << table_header << '\n';
    bench_range<MIN_HEIGHT, MAX_HEIGHT, Hash>(name);
    log_file << '\n';
}

int main()
{
    fs::create_directories("./log");
    std::string timestamp = std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(
                                               std::chrono::system_clock::now().time_since_epoch())
                                               .count());
    std::string log_file_name = std::string("./log/benchmark_native_") + timestamp +
                                std::string(".log");
    log_file.open(log_file_name);

    std::cout << "Logging to " << log_file_name << "...\n";

    log_file << std::boolalpha;
    libff::inhibit_profiling_info = true;
    libff::inhibit_profiling_counters = true;

    libff::default_ec_pp::init_public_params();

    log_file << "Hash Functions Benchmark"
             << "\n";
    log_file << "Prime:\t" << FieldT::mod << "\n";

#ifdef MULTICORE
    omp_set_num_threads(NUM_THREADS);
    log_file << "Threads:\t" << omp_get_max_threads() << "\n";
#else
    log_file << "Threads:\t1\n";
#endif

    bench<Sha256>("SHA-256 (2:1)");

    bench<Arion<FieldT, 2, 1, 6>>("Arion");
    bench<Arion<FieldT, 4, 1, 5>>("Arion");
    bench<Arion<FieldT, 8, 1, 4>>("Arion");

    bench<Anemoi<FieldT, 2, 2, 14>>("Anemoi");
    bench<Anemoi<FieldT, 4, 2, 12>>("Anemoi");
    bench<Anemoi<FieldT, 8, 2, 11>>("Anemoi");

    bench<Griffin<FieldT, 2, 1, 12>>("Griffin");
    bench<Griffin<FieldT, 4, 4, 9>>("Griffin");
    bench<Griffin<FieldT, 8, 4, 9>>("Griffin");

    bench<Poseidon<FieldT, 2, 1, 4, 55>>("Poseidon");
    bench<Poseidon<FieldT, 4, 1, 4, 56>>("Poseidon");
    bench<Poseidon<FieldT, 8, 1, 4, 56>>("Poseidon");

    bench<Poseidon2<FieldT, 2, 4, 55>>("Poseidon2");
    bench<Poseidon2<FieldT, 4, 4, 56>>("Poseidon2");
    bench<Poseidon2<FieldT, 8, 4, 56>>("Poseidon2");

    bench<Rescue<FieldT, 2, 1, 14>>("Rescue");
    bench<Rescue<FieldT, 4, 1, 9>>("Rescue");
    bench<Rescue<FieldT, 8, 1, 8>>("Rescue");

    log_file.close();

    return 0;
}
