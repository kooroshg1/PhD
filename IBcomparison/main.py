__author__ = 'koorosh'
import numpy as np
import matplotlib.pyplot as plt

def solve_poisson(x, n, BC):
    def generate_stiffness(x, n):
        A = np.zeros((n+2, n+2))
        for ni in range(1, n+1):
            dxl = x[ni] - x[ni-1]
            dxr = x[ni+1] - x[ni]
            A[ni, ni-1] = dxr
            A[ni, ni] = -(dxl + dxr)
            A[ni, ni+1] = dxl
            A[ni, :] /= (dxr * dxl * (dxr + dxl) * 0.5)
        return A


    def generate_force(A, BC, n):
        bcl, bcr = BC
        F = np.zeros((n+2, 1))
        F[1, 0] = -A[1, 0] * bcl
        F[n, 0] = -A[n, n+1] * bcr
        return F


    K = generate_stiffness(x, n)
    F = generate_force(K, BC, n)

    K = K[1:n+1, 1:n+1]
    F = F[1:n+1]

    x = np.linalg.solve(K, F)
    x = np.append(BC[0], x)
    x = np.append(x, BC[1])
    return x

n = 3

x1 = np.linspace(0.0, 1.0, n + 2)
x2 = np.linspace(0.0, 1.0, n + 2)
x2[-1] = x2[-1] - 0.1
x3 = np.linspace(0.0, 0.9, n + 2)


y1 = solve_poisson(x1, n, [0.0, 1.0])
y2 = solve_poisson(x2, n, [0.0, 1.0])
y3 = solve_poisson(x3, n, [0.0, 1.0])

print(y1)
print(y2)
plt.figure()
plt.plot(x1, y1, 'k-+',
         x2, y2, 'wo',
         x3, y3, 'r-+',
         ms=10, mew=2.0)
plt.show()