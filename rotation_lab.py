"""
SymPy로 회전 행렬·쿼터니언·exp(A+B) vs R_B R_A 수치 확인용.
실행: .venv/bin/python rotation_lab.py
"""
from sympy import Matrix, pi, cos, sin, sqrt, simplify, N
import sympy as sp


def Rx(a):
    return Matrix([[1, 0, 0], [0, cos(a), -sin(a)], [0, sin(a), cos(a)]])


def Ry(a):
    return Matrix([[cos(a), 0, sin(a)], [0, 1, 0], [-sin(a), 0, cos(a)]])


def hat(v):
    vx, vy, vz = v[0], v[1], v[2]
    return Matrix([[0, -vz, vy], [vz, 0, -vx], [-vy, vx, 0]])


def rot_from_axis_angle(axis, angle):
    ax = axis / sqrt(axis.dot(axis))
    return Matrix.exp(angle * hat(ax))


def frob(M):
    return sqrt((M.T * M).trace())


def main():
    deg90 = pi / 2
    deg45 = pi / 4

    # --- 1. 비가환: Rx(90°) Ry(90°) vs Ry(90°) Rx(90°) ---
    Rxy = simplify(Rx(deg90) * Ry(deg90))
    Ryx = simplify(Ry(deg90) * Rx(deg90))
    print("=== 1. Rx*Ry vs Ry*Rx (90°) ===")
    print("Rx*Ry:\n", Rxy)
    print("Ry*Rx:\n", Ryx)
    same = [[bool(Rxy[i, j] == Ryx[i, j]) for j in range(3)] for i in range(3)]
    print("원소 일치 여부 (행렬):", same)

    # --- 2. z축 45° 쿼터니언 (w, x, y, z) ---
    th = deg45
    w, x, y, z = cos(th / 2), 0, 0, sin(th / 2)
    print("\n=== 2. z축 45° (w,x,y,z) ===")
    print("w = cos(θ/2) ≈", N(w, 10))
    print("(x,y,z) = sin(θ/2)*(0,0,1) → z ≈", N(z, 10))

    # --- 3. exp(A+B) vs R_B R_A (A=θ hat(ez), B=θ hat(ex)) ---
    ez = Matrix([0, 0, 1])
    ex = Matrix([1, 0, 0])

    for label, theta in [("0.1 rad", sp.Rational(1, 10)), ("1.5 rad", sp.Rational(3, 2))]:
        A = theta * hat(ez)
        B = theta * hat(ex)
        R_A = rot_from_axis_angle(ez, theta)
        R_B = rot_from_axis_angle(ex, theta)
        prod = simplify(R_B * R_A)
        ssum = simplify(Matrix.exp(A + B))
        diff = prod - ssum
        comm = simplify(A * B - B * A)
        print(f"\n=== 3. {label} (R_B R_A vs exp(A+B)) ===")
        print("||R_B R_A - exp(A+B)||_F ≈", N(frob(diff), 12))
        print("||(1/2)[A,B]||_F ≈", N(frob(comm / 2), 12))


if __name__ == "__main__":
    main()
