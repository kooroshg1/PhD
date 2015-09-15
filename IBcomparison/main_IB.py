__author__ = 'koorosh'
import numpy as np
import matplotlib.pyplot as plt

n = 3
x = np.linspace(0.0, 1.0, n + 2)

xBC = 0.75 + 0.25 / 2.0
nBC = x - xBC
m = min(i for i in nBC if i >= 0)
if m == 0:
    isOnNode = True
else:
    isOnNode = False
endIndex = np.where(nBC == m)[0][0]

# Calculating the weights
if ~isOnNode:
    leftWeight = 1 - (xBC - x[endIndex - 1]) / (x[endIndex] - x[endIndex - 1])
    rightWeight = (xBC - x[endIndex - 1]) / (x[endIndex] - x[endIndex - 1])
else:
    leftWeight, rightWeight = 0.0, 0.0

def laplace_operator(x, n):
    A = np.zeros((n+2, n+2))
    dx = x[2] - x[1]
    for ni in range(1, n+1):
        A[ni, ni-1] = 1
        A[ni, ni] = -2
        A[ni, ni+1] = 1
    return A


def discrete_force(endIndex, n, weight):
    if len(weight) == 2:
        leftWeight, rightWeight = weight
    else:
        weight = weight
    Af = np.zeros((n+2, n+2))
    Af[endIndex-1, endIndex-1] = leftWeight
    Af[endIndex-1, endIndex] = rightWeight
    return Af

def boundary_condition(A, n, BC):
    bcl, bcr = BC
    F = np.zeros((n+2, 1))
    F[1, 0] = -A[1, 0] * bcl
    F[n, 0] = -A[n, n+1] * bcr
    return F

A = laplace_operator(x, n) + discrete_force(endIndex, n, [leftWeight, rightWeight])
b = boundary_condition(A, n, [0.0, 1.0])
print(laplace_operator(x, n))
print()
print(discrete_force(endIndex, n, [leftWeight, rightWeight]))
print()
print(boundary_condition(A, n, [0.0, 1.0]))
print()

A = A[1:n+1, 1:n+1]
b = b[1:n+1]
sol = np.linalg.solve(A, b)
print(sol)
'''
A = np.array([[0.0, 0.0, 0.0, 0.0, 0.0],
              [1.0, -2.0, 1.0, 0.0, 0.0],
              [0.0, 1.0, -2.0, 1.0, 0.0],
              [0.0, 0.0, 1.0, -3.0, 0.0],
              [0.0, 0.0, 0.0, 0.0, 0.0]])


F = np.matrix([0.0, 0.0, 0.0, -2.0, 0.0])
F = F.transpose()

A = A[1:4, 1:4]
F = F[1:4]

sol = np.linalg.solve(A, F)
sol = np.append(0, sol)
sol = np.append(sol, 1.0)

plt.figure()
plt.plot([0.0, 0.25, 0.5, 0.75, 0.875], sol)
plt.show()
'''