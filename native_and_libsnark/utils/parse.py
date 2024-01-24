from sys import argv
from math import log2


def f_calls(t, h):
    return (t**h - 1) // (t - 1)


def main():
    if len(argv) < 2:
        print(f"Usage: {argv[0]} <log_path>")
        exit(-1)

    fname = argv[1]
    with open(fname, "rt") as log:
        lines = log.readlines()
        if fname.find("_ops_") != -1:
            for i in range(len(lines)):
                if lines[i].startswith("Operation"):
                    print(lines[i - 1], end="")
                    for j in range(i + 3, i + 13):
                        kops = float(lines[j].split("\t")[2])
                        ns_per_op = 1_000_000 / kops
                        print(f"{lines[j].split()[0]}, {ns_per_op:.2f}")
                    print()
        elif fname.find("_native_") != -1 and False:
            time = []
            for i in range(len(lines)):
                if lines[i].startswith("Height"):
                    name = lines[i - 1].strip()
                    data = list(map(float, lines[i + 3].strip().split("\t")))
                    print(f"{name}: ", end="")
                    if data[0] == 18:
                        c = f_calls(2, data[0])
                    elif data[0] == 9:
                        c = f_calls(4, data[0])
                    elif data[0] == 6:
                        c = f_calls(8, data[0])
                    else:
                        c = 1
                    time.append(data[1] / c * 1000)
                    if name.find("-DM") != -1:
                        if name.find("110") != -1:
                            speedup = time[-4] / time[-1]
                        else:
                            speedup = time[-5] / time[-1]
                        print(f"{time[-1]:.2f} ({speedup:.2f}x)")
                    else:
                        print(f"{time[-1]:.2f}")
        elif fname.find("_native_") != -1:
            for i in range(len(lines)):
                if lines[i].startswith("Height"):
                    name = lines[i - 1].strip()
                    print(name)
                    j = i + 1
                    base_data = 1
                    while lines[j][0].isnumeric():
                        data = lines[j].split("\t")
                        data = [int(data[0]), float(data[1])]
                        if j == i + 1:
                            base_data = data[1]
                        data[0] = 4 ** data[0]
                        data[1] = data[1]
                        print(f"({data[0]},{data[1]})")
                        j += 1


if __name__ == "__main__":
    main()
