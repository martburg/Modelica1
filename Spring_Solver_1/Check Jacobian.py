import numpy as np
import os

N = 57
THRESH = 1e-6

def read_jacobian_from_file(filename):
    data = []
    with open(filename, 'r') as f:
        for line in f:
            parts = line.strip().split()
            data.extend([float(x) for x in parts])
    if len(data) != N * N:
        raise ValueError(f"Expected {N*N} elements, got {len(data)}")
    return np.array(data).reshape((N, N))

def check_singular_rows(J):
    for i in range(N):
        row_sum = np.sum(np.abs(J[i, :]))
        if row_sum < THRESH:
            print(f"⚠️  Row {i} is nearly zero! sum = {row_sum:.3e}")

if __name__ == "__main__":
    print("Current working directory:", os.getcwd())
    assert os.path.exists(".\Spring_Solver_1\jacobian.txt"), "❌ File not found: jacobian.txt"
    filename = ".\Spring_Solver_1\jacobian.txt"  # Replace with your filename
    J = read_jacobian_from_file(filename)
    check_singular_rows(J)