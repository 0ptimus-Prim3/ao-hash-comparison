#include "hash/poseidon/poseidon.hpp"
#include "util/measure.hpp"
#include "util/string_utils.hpp"
#include <cstring>
#include <iostream>

using Hash = Poseidon<libff::Fr<libff::default_ec_pp>, 2, 1>;

static bool run_tests()
{
    auto msg = HEXARRAY(
        00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001);
    uint8_t dig[Hash::DIGEST_SIZE]{};
    auto real_dig = HEXARRAY(d255a01dc046c46df0ef710883852fb60ebcf1749480a17789d67a42311c7e66);

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
    std::cout << "\n==== Testing Poseidon ====\n";

    bool all_check = run_tests();

    std::cout << "\n==== " << (all_check ? "ALL TESTS SUCCEEDED" : "SOME TESTS FAILED")
              << " ====\n\n";

#ifdef MEASURE_PERFORMANCE
#endif

    return 0;
}
