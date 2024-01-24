from math import floor, ceil, log2, comb as binomial
from sys import argv

BLS381 = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
BN254 = 0x30644e72e131a029b85045b68181585d2833e84879b9709143e1f593f0000001

def get_number_of_rounds(t, d, security_level, p=BN254):
    sl = 1 << (security_level >> 1)
    nr = ceil((2.5 * security_level) / (log2(p) - log2(d - 1)))
    rgb = 0
    while rgb < 100:
        c = min(binomial(rgb * (d + t) + 1, 1 + t * rgb), binomial(d ** rgb + 1 + rgb, 1 + rgb))
        if c >= sl:
            break
        rgb += 1

    return ceil(1.2 * max(6, nr, 1 + rgb))


def get_r1cs_size(t, d, nr):
    addchain = (0, 0, 1, 2, 2, 3, 3, 4, 3, 4)

    return (2 * addchain[d] + (t - 2) * 2) * nr


def main():
    if len(argv) < 4:
        print(f"Usage: python {argv[0]} branch d security_level")
        exit(1)

    t = int(argv[1])
    d = int(argv[2])
    security_level = int(argv[3])
    nr = get_number_of_rounds(t, d, security_level)
    print(f"Number of rounds: {nr}")
    r1cs = get_r1cs_size(t, d, nr)
    print(f"R1CS size: {r1cs}")


if __name__ == '__main__':
    main()
