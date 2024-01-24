from math import floor, ceil, comb as binomial
from sys import argv


def get_number_of_rounds(m, capacity, alpha, security_level):
    # get number of rounds for Groebner basis attack
    rate = m - capacity
    def dcon(N): return floor(0.5 * (alpha - 1) * m * (N - 1) + 2)
    def v(N): return m * (N - 1) + rate
    target = 2 ** security_level
    for l1 in range(1, 25):
        if binomial(v(l1) + dcon(l1), v(l1)) ** 2 > target:
            break
    # set a minimum value for sanity and add 50%
    return ceil(1.5 * max(5, l1))


def get_r1cs_size(t, d, nr):
    addchain = (0, 0, 1, 2, 2, 3, 3, 4, 3, 4)

    return 2 * addchain[d] * t * nr


def main():
    if len(argv) < 5:
        print(f"Usage: python {argv[0]} rate capacity alpha security_level")
        exit(1)

    rate = int(argv[1])
    capacity = int(argv[2])
    alpha = int(argv[3])
    security_level = int(argv[4])
    m = rate + capacity

    nr = get_number_of_rounds(m, capacity, alpha, security_level)
    print(f"Number of rounds: {nr}")

    r1cs = get_r1cs_size(m, alpha, nr)
    print(f"R1CS size: {r1cs}")


if __name__ == '__main__':
    main()
