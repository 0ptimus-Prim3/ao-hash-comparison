from curses.ascii import isdigit
from matplotlib import pyplot as plt
import numpy as np
import matplotlib
from sys import argv


def main():
    if len(argv) < 2:
        print("Usage: plot.py <log_file>")
        exit(1)

    with open(argv[1], "r") as f:
        lines = f.readlines()
        i = 0
        fig_i = 0
        colors = ["r", "b", "g", "y", "m", "c"]
        color = 0
        old_title = None

        matplotlib.use("QtAgg")

        plt.ion()
        plt.figure(argv[1])

        ax = plt.subplot()
        ax.xaxis.set_major_locator(plt.MultipleLocator(4))
        ax.yaxis.set_major_locator(plt.MultipleLocator(1))
        data = []
        while i < len(lines):
            if lines[i].startswith("/*"):
                while i < len(lines) and not lines[i].startswith("*/"):
                    i += 1
            if lines[i].startswith("Height"):
                data.append(np.empty(0))
                title = lines[i-1]
                metrics = lines[i].split()
                i += 1
                while i < len(lines) and isdigit(lines[i][0]):
                    data[-1] = np.append(data[-1], np.asarray(lines[i].split()).astype(float))
                    i += 1
                data[-1] = data[-1].reshape(-1, len(metrics)).T
                col = 1
                if old_title:
                    color += not title.startswith(old_title.split()[0])
                old_title = title

                for j in range(col, col+1):
                    ax.plot(data[-1][0, :], np.log(data[-1][j, :]),
                             label=f"{title}", marker="o")
            else:
                i += 1
        plt.legend()
        plt.xlabel("# Cores")
        plt.ylabel("Time (s)")
        plt.grid()
        plt.ioff()
        plt.show()


if __name__ == "__main__":
    main()
