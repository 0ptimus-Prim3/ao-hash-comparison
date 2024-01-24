from math import floor, ceil, comb
from sys import argv


def get_number_of_rounds_old2(ell, alpha, security_level):
    sl = 1 << security_level
    for r in range(100):
        b = comb(2 * ell * r + alpha + 1 + 2 * (ell * r - 2), 2 * ell * r)
        if b * b >= sl:
            return 1 + ell + r

    return 10


def get_number_of_rounds_old1(ell, alpha, security_level):
    kappa = (1, 1, 1, 1, 2, 2, 4, 4, 7, 7)
    sl = 1 << security_level

    for r in range(100):
        k = kappa[alpha] if alpha < len(kappa) else 9
        c_alg = comb(4 * ell * r + k, 2 * ell * r)
        c_alg *= c_alg
        if c_alg >= sl:
            return max(8, min(5, 1 + ell) + 2 + r)

    return 8


def get_number_of_rounds(ell, alpha, security_level):
    kappa = (1, 1, 1, 1, 2, 2, 4, 4, 7, 7)
    sl = 1 << security_level

    for r in range(100):
        k = kappa[alpha] if alpha < len(kappa) else 9
        c_alg = comb(2 * ell * r + alpha + 1 + 2 * (ell * r - 2), 2 * ell * r)
        c_alg *= c_alg
        if c_alg >= sl:
            return max(10, r + 1 + ell)

    return 10


def get_r1cs_size(ell: int, alpha: int, nr: int) -> int:
    addchain = (0, 0, 1, 2, 2, 3, 3, 4, 3, 4)

    if alpha >= len(addchain):
        return -1

    return (addchain[alpha] + 2) * ell * nr


def main():
    if len(argv) < 4:
        print(f"Usage: python {argv[0]} l alpha security_level")
        exit(1)

    ell = int(argv[1])
    alpha = int(argv[2])
    security_level = int(argv[3])
    nr = get_number_of_rounds(ell, alpha, security_level)
    print(f"Number of rounds: {nr}")
    r1cs = get_r1cs_size(ell, alpha, nr)
    print(f"R1CS size: {r1cs}")


if __name__ == '__main__':
    main()
