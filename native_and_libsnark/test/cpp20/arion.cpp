#include "hash/arion/arion.hpp"
#include "util/string_utils.hpp"
#include <cstring>
#include <iostream>

using ppT = libff::default_ec_pp;
using FieldT = libff::Fr<ppT>;
using Hash = Arion<FieldT, 2, 1>;

static bool run_tests()
{
    //There are no test vectors for MiMC, so we assume our implementation to be correct
    uint8_t msg[Hash::BLOCK_SIZE]{};
    uint8_t dig[Hash::DIGEST_SIZE]{};

#ifdef CURVE_ALT_BN128
    auto real_dig = HEXARRAY(283d55b60454ebe917f242ba13d0256c0df7d4f52659d287b04adc467e28c389);
#else
    auto real_dig = HEXARRAY(67d43da37c30152ff5aacb2d7d4940748904abe5dfc77f637109a2a96b178de1);
#endif

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
    std::cout << "\n==== Testing Arion ====\n";

    bool all_check = run_tests();

    std::cout << "\n==== " << (all_check ? "ALL TESTS SUCCEEDED" : "SOME TESTS FAILED")
              << " ====\n\n";

#ifdef MEASURE_PERFORMANCE
#endif

    return 0;
}
