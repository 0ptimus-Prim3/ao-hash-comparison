from math import ceil, log
from sys import argv


def get_number_of_rounds(p, d):
    return ceil(log(p, d))


def get_r1cs_size(d, n):
    addchain = (0, 0, 1, 2, 2, 3, 3, 4, 3, 4)

    return addchain[d] * n


def main():
    if len(argv) < 3:
        print(f"Usage: python {argv[0]} p d")
        exit(1)

    p = int(argv[1])
    d = int(argv[2])
    n = get_number_of_rounds(p, d)
    print(f"Number of rounds: {n}")
    r1cs = get_r1cs_size(d, n)
    print(f"R1CS size: {r1cs}")


if __name__ == '__main__':
    main()
