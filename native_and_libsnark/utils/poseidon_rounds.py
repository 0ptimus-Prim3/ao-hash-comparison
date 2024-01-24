from math import ceil, log2, comb
from sys import argv

BLS381 = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
BN254 = 0x30644e72e131a029b85045b68181585d2833e84879b9709143e1f593f0000001


def get_number_of_rounds(t, alpha, security_level):
    base = {80: 35, 128: 56, 256: 111}

    br = base[security_level]

    if security_level == 128:
        if t == 3:
            return 8, 57
        if t == 5:
            return 8, 60
        if t == 9:
            return 8, 63

    l_alpha = ceil(log2(t)/log2(alpha))

    return 8, ceil(1.075 * (l_alpha + br))


def get_number_of_rounds_new(t, alpha, security_level, p=BLS381):
    sl = min(security_level, log2(p))
    min_sec = sl/log2(alpha) + log2(t)/log2(alpha) - 5

    grob_sec_1 = 0
    qty = 0
    limit = 2 ** (security_level // 2)

    # In Sponge mode we assume r = t - 1, c = 1
    while qty < limit:
        grob_sec_1 += 1
        upper = 4 * t + 5 * (t - 1) + 2 * grob_sec_1 + alpha
        lower = 3 * (t - 1) + grob_sec_1 + alpha
        qty = comb(upper, lower)


    grob_sec_2 = 1/log2(alpha) * sl - 6
    grob_sec_3 = t - 7 + 1/log2(alpha) * \
        min(security_level / (t + 1), log2(p)/2)
    grob_sec_4 = security_level / (2 * log2(alpha)) - 5 * t + 4
    grob_sec = max(grob_sec_1, grob_sec_2, grob_sec_3, grob_sec_4)

    return 8, ceil(1.075 * ceil(max(min_sec, grob_sec)))


def get_r1cs_size(t, d, nf, np):
    addchain = (0, 0, 1, 2, 2, 3, 3, 4, 3, 4)

    return addchain[d] * (nf + np) + addchain[d] * (t - 1) * nf


def main():
    if len(argv) < 4:
        print(f"Usage: python {argv[0]} branches alpha security_level")
        exit(1)

    t = int(argv[1])
    alpha = int(argv[2])
    security_level = int(argv[3])
    nf, np = get_number_of_rounds(t, alpha, security_level)
    print(f"Number of rounds: {nf, np}")
    nf, np = get_number_of_rounds_new(t, alpha, security_level)
    print(f"Number of rounds new: {nf, np}")
    r1cs = get_r1cs_size(t, alpha, nf, np)
    print(f"R1CS size: {r1cs}")


if __name__ == '__main__':
    main()
