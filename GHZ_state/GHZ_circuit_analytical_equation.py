#!/usr/bin python3
'''
the GHZ circuit equation - as calculated by hand with input qubits 0,0,0
WARNING: This is wrong!
(does not agree with first column of GHZ_circuit_unitary nor fig3_rho from GHZ_circuit_with_qutip
'''

from numpy import cos, sin, pi, random, around
from pprint import pformat

i = 1j

phis = random.uniform(0, 2*pi, 6)
phis = [pi/3] * 6
phi_1, phi_2, phi_3, phi_4, phi_5, phi_6 = phis


a = phi_1 / 2
b = phi_2 / 2
c = phi_3 / 2

# -----------------------------------------------------
# The 8x8 state vector after the 3 Rx's and the 2 CNOTs
# -----------------------------------------------------

# NOTE: here the 0s correspond to cosines, and the 1s to sines!
# CNOT(a,b) causes the qubit b to flip if and only if the qubit a is 1
# this half was untouched by the CNOTs as they are rooted on the first qubit (which is always zero here)
S = +1 * cos(a) * cos(b) * cos(c)  # 0 = 000
T = -i * cos(a) * cos(b) * sin(c)  # 1 = 001
U = -i * cos(a) * sin(b) * cos(c)  # 2 = 010
V = -1 * cos(a) * sin(b) * sin(c)  # 3 = 011

# This half was reversed order due to being passed through -> CNOT(0,1) -> CNOT(0,2)
# this half is identical to the first half if sin <-> cosine, and multiplied by i
# this means a phase shift of 90º and a rotation by 90º anticlockwise in the imaginary plane
W = +i * sin(a) * sin(b) * sin(c)  # 4 = 100 -> 110 -> 111
X = -1 * sin(a) * sin(b) * cos(c)  # 5 = 101 -> 111 -> 110
Y = -1 * sin(a) * cos(b) * sin(c)  # 6 = 110 -> 100 -> 101
Z = -i * sin(a) * cos(b) * cos(c)  # 7 = 111 -> 101 -> 100

# --------------------------------------------------
# Derived from Rx(phi_4) (x) Rx(phi_5) (x) Rx(phi_6)
# --------------------------------------------------

# these angles combine when finding the tensor product of the Rx(phi_4) and Rx(phi_5) gates
# (after the entangling of the 3 qubits by the CNOT gates)
# due to the classic trigonometric rules:
# cosAcosB - sinAsinB = cos(A+B)
# sinAcosB + cosAsinB = sin(A+B)
# proved using Euler's formula: https://math.stackexchange.com/a/860634
d = (phi_4 + phi_5) / 2
e = phi_6 / 2

# Derived from Rx(phi_4) (x) Rx(phi_5) (x) Rx(phi_6)
# ((A, B),(C, D)) is the 2x2 structure repeated on the diagonal of the 8x8 operation
A = +1 * cos(d) * cos(e) +i * sin(d) * sin(e)
B = -1 * cos(d) * sin(e) -i * sin(d) * cos(e)
C = +1 * cos(d) * sin(e) -i * sin(d) * cos(e)
D = +1 * cos(d) * cos(e) -i * sin(d) * sin(e)

# Applying the derived Rx,Rx,Ry 8x8 matrix to the Rx,Rx,Rx + CNOT'd 8 state column vector
GHZ = (
    A*S + B*T,
    C*S + D*T,
    A*U + B*V,
    C*U + B*V,
    A*W + B*X,
    C*W + D*X,
    A*Y + B*Z,
    C*Y + D*Z,
)
GHZ = [around(state, 3) for state in GHZ]

if __name__ == '__main__':
    print('phis (in degrees) = ', [round(phi * (180 / pi), 1) for phi in phis])
    print('GHZ rho for |000> input:\n', pformat(GHZ))
