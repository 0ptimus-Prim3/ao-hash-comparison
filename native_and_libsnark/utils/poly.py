import numpy as np

n = 3
N = 1 << n
z1 = np.zeros((N, N),dtype=np.int64)
z2 = np.zeros((N, N),dtype=np.int64)

for i in range(N):
    for j in range(N):
        z1[i,j] = i ^ j
        if i % 2:
            z2[i,j] = abs(i - j) % N
        else:
            z2[i,j] = (i + j) % N

print(z1)
print()
print("=================================")
print()
print(z2)
