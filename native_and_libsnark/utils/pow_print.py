from sys import argv

# compute fast exponentiation non recursively


def fast_pow(n, exp) -> int:
    result = 1
    while (exp > 0):
        if (exp & 1):
            result *= n
        n *= n
        exp >>= 1
    return result


def fast_pow_print_steps(exp):
    i = 1
    while exp > 0:
        if exp % 2 == 1:
            if i != 1:
                print(f"x *= t{';' if i % 10 == 0 else ','}")
                i += 1
        exp //= 2
        if exp > 0:
            print(f"t *= t{';' if i % 10 == 0 else ','}")
            i += 1


def main():
    if (len(argv) < 2):
        print("Usage: pow_min.py <exp>")
        exit(1)

    fast_pow_print_steps(int(argv[1], base=0))


if __name__ == "__main__":
    main()
