#include "hash/poseidon2/poseidon2.hpp"
#include "util/measure.hpp"
#include "util/string_utils.hpp"
#include <cstring>
#include <iostream>

using Hash = Poseidon2<libff::Fr<libff::default_ec_pp>, 2>;

static bool run_tests()
{
    auto msg = HEXARRAY(
        00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000);
    uint8_t dig[Hash::DIGEST_SIZE]{};
    auto real_dig = HEXARRAY(e58ee5fb8b73692328321ff1bfb90fb4e19833220d41a9dbb3ae7263b192d908);

    bool check = true;
    bool all_check = true;

    std::cout << std::boolalpha;

    std::cout << "Hashing libff... ";
    check = true;
    Hash::hash_oneblock(dig, msg.data());
    check = memcmp(dig, real_dig.data(), sizeof(dig)) == 0;
    std::cout << check << '\n';
    all_check &= check;

    return all_check;
}

int main()
{
    std::cout << "\n==== Testing Poseidon2 ====\n";

    bool all_check = run_tests();

    std::cout << "\n==== " << (all_check ? "ALL TESTS SUCCEEDED" : "SOME TESTS FAILED")
              << " ====\n\n";

#ifdef MEASURE_PERFORMANCE
#endif

    return 0;
}
