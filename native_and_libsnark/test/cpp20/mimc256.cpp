#include "hash/mimc/mimc256.hpp"
#include "util/string_utils.hpp"
#include <cstring>
#include <iostream>

using ppT = libff::default_ec_pp;
using FieldT = libff::Fr<ppT>;
using Hash = Mimc256<FieldT>;

static bool run_tests()
{
    //There are no test vectors for MiMC, so we assume our implementation to be correct
    uint8_t msg[Hash::BLOCK_SIZE]{};
    uint8_t dig[Hash::DIGEST_SIZE]{};
    auto real_dig = HEXARRAY(0ce0fabfe10dffc378fbf7d413c1cd6c63d2114abcef1dfc6680223be4523321);

    bool check = true;
    bool all_check = true;

    std::cout << std::boolalpha;

    std::cout << "Hashing... ";
    check = true;

    Hash::hash_oneblock(dig, msg);
    check = memcmp(dig, real_dig.data(), sizeof(dig)) == 0;

    std::cout << hexdump(dig) << '\n';
    std::cout << check << '\n';
    all_check &= check;

    return all_check;

   return true;
}

int main()
{
    std::cout << "\n==== Testing Hash ====\n";

    bool all_check = run_tests();

    std::cout << "\n==== " << (all_check ? "ALL TESTS SUCCEEDED" : "SOME TESTS FAILED")
              << " ====\n\n";

#ifdef MEASURE_PERFORMANCE
#endif

    return 0;
}
