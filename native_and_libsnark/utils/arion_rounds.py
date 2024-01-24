from math import log2, sqrt, gcd
from sys import argv

BLS381 = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
BN254 = 0x30644e72e131a029b85045b68181585d2833e84879b9709143e1f593f0000001


def get_number_of_rounds(n, d1, d2, security_level, p=BN254):
    w = 2
    r = 4
    lp = log2(p)

    while r < 1000:
        d = (2 ** (n - 1) * d2 * (d1 + 1) - d1 * d2) ** r
        ld = log2(d)
        lld = log2(ld)
        c_common = (d * lp * ld * lld) + (d * ld * ld * lld)
        c_det = (n * d ** w) + c_common
        c_prob = (sqrt(n) * d ** (2 + (n - 1)/n)) + c_common
        c = min(c_det, c_prob)
        lc = log2(c)
        if round(lc) >= security_level:
            break
        r += 1

    return lc, r


def get_r1cs_size(t, d1, d2: int, nr):
    addchain = (0, 0, 1, 2, 2, 3, 3, 4, 3, 4)

    # we assume d2 = 2^n + 1
    return (d2.bit_length() + (t - 1) * (addchain[d1] + 2)) * nr


def main():
    if len(argv) < 5:
        print(f"Usage: python {argv[0]} branch d1 d2 security_level")
        exit(1)

    t = int(argv[1])
    d1 = int(argv[2])
    d2 = int(argv[3])
    security_level = int(argv[4])

    if gcd(BLS381 - 1, d2) != 1:
        print("Warning: x^d2 is not a permutation in BLS381")
    if gcd(BN254 - 1, d2) != 1:
        print("Warning: x^d2 is not a permutation in BN254")

    lc, nr = get_number_of_rounds(t, d1, d2, security_level)
    print(f"Number of rounds: {nr} (security level: {round(lc, 2)})")
    r1cs = get_r1cs_size(t, d1, d2, nr)
    print(f"R1CS size: {r1cs}")


if __name__ == '__main__':
    main()
