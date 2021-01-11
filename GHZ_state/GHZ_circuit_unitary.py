#!/usr/bin python3
'''
GHZ circuit equation (general input, not just |000>)
'''

from numpy import cos, sin, pi, random, array, eye, kron, trace, around, sqrt, transpose
from pprint import pformat


i = 1j

# phis = random.uniform(0, 2*pi, 6)
# phis = [pi/3] * 6
# phis = [pi, pi/2, pi/3, pi/4, pi/5, pi/6]
phis_deg = [69.1, 220.1, 131.7, 160.3, 188.8, 80.9]
phis = [phi * (pi/180) for phi in phis_deg]
phi_1, phi_2, phi_3, phi_4, phi_5, phi_6 = phis

# 2x2 identity matrix
I = eye(2)
I8 = eye(8)

# Pauli matrices
sx = array([[0, 1], [1, 0]])
sy = array([[0, -1j], [1j, 0]])
sz = array([[1, 0], [0, -1]])

zero = array([[1],[0]])  # aka 'spin up' or |0>
one = array([[0],[1]])  # aka 'spin down' or |1>


def triple_kron(matrix_1, matrix_2, matrix_3):
    ''' Does m1 (x) m2 (x) m3 '''
    return kron(kron(matrix_1, matrix_2), matrix_3)


triple_0 = triple_kron(zero, zero, zero)  # |000>
triple_1 = triple_kron(one, one, one)  # |111>
GHZ_ket = (1 / sqrt(2)) * (triple_0 + triple_1)
GHZ_bra = GHZ_ket.conj().T


def Rx(phi):
    a = phi / 2
    return array([
        [cos(a), -i * sin(a)],
        [-i * sin(a), cos(a)],
    ])


def Ry(phi):
    a = phi / 2
    return array([
        [cos(a), -sin(a)],
        [sin(a), cos(a)],
    ])


# if CNOT_AB qubit A has value 1, flip qubit 2
CNOT_01 = triple_kron(zero @ zero.conj().T, I, I) + triple_kron(one @ one.conj().T, sx, I)
CNOT_02 = triple_kron(zero @ zero.conj().T, I, I) + triple_kron(one @ one.conj().T, I, sx)


Rx_1 = triple_kron(Rx(phi_1), I, I)
Rx_2 = triple_kron(I, Rx(phi_2), I)
Rx_3 = triple_kron(I, I, Rx(phi_3))
Rx_4 = triple_kron(Rx(phi_4), I, I)
Rx_5 = triple_kron(I, Rx(phi_5), I)
Ry_6 = triple_kron(I, I, Ry(phi_6))


GHZ_circuit_unitary = Ry_6 @ Rx_5 @ Rx_4 @ CNOT_02 @ CNOT_01 @ Rx_3 @ Rx_2 @ Rx_1
GHZ_circuit_unitary_rounded = around(GHZ_circuit_unitary, 3)

# NOTE: this transposes the bottom right 4x4 of the 8x8 unitary
# GHZ_circuit_unitary_reverse = Rx_1 @ Rx_2 @ Rx_3 @ CNOT_01 @ CNOT_02 @ Rx_4 @ Rx_5 @ Ry_6
# GHZ_circuit_unitary_reverse_rounded = around(GHZ_circuit_unitary_reverse, 3)

Fidelity_bra_rho_ket = GHZ_bra @ GHZ_circuit_unitary @ GHZ_ket
Fidelty_wiki = (GHZ_bra @ GHZ_circuit_unitary[:,0])
print('huh', Fidelty_wiki, abs(Fidelty_wiki), abs(Fidelty_wiki) ** 2)
Fidelty_wiki = abs(GHZ_bra @ GHZ_circuit_unitary) ** 2

print('huh', around(GHZ_circuit_unitary.conj().T, 2))
print('huh2', around(GHZ_circuit_unitary @ GHZ_circuit_unitary.conj().T, 2))
rho = GHZ_circuit_unitary @ GHZ_ket @ GHZ_bra @ GHZ_circuit_unitary.conj().T
rho = GHZ_circuit_unitary @ triple_0 @ triple_0.conj().T @ GHZ_circuit_unitary.conj().T
print('rho', rho)
Fidelity_bra_rho_ket = GHZ_bra @ rho @ GHZ_ket
Fidelty_wiki = abs(GHZ_bra @ rho) ** 2


if __name__ == '__main__':
    print('phis (in degrees) = ', [round(phi * (180 / pi), 1) for phi in phis])
    # print('GHZ unitary:\n', pformat(GHZ_circuit_unitary_rounded))
    # big endian
    print('GHZ |000> column:\n', pformat(GHZ_circuit_unitary_rounded[:,0]))
    print('GHZ |001> column:\n', pformat(GHZ_circuit_unitary_rounded[:,1]))
    print('GHZ |010> column:\n', pformat(GHZ_circuit_unitary_rounded[:,2]))
    print('GHZ |011> column:\n', pformat(GHZ_circuit_unitary_rounded[:,3]))
    print('GHZ |100> column:\n', pformat(GHZ_circuit_unitary_rounded[:,-4]))
    print('GHZ |101> column:\n', pformat(GHZ_circuit_unitary_rounded[:,-3]))
    print('GHZ |110> column:\n', pformat(GHZ_circuit_unitary_rounded[:,-2]))
    print('GHZ |111> column:\n', pformat(GHZ_circuit_unitary_rounded[:,-1]))

    # print('GHZ |000> column (reverse):\n', pformat(GHZ_circuit_unitary_reverse_rounded[:,0]))
    print('Fidelity <GHZ|rho|GHZ> = ', pformat(around(Fidelity_bra_rho_ket, 3)))
    print('Fidelity |<GHZ|rho_000>|^2 = ', pformat(around(Fidelty_wiki, 3)))
    print('Fidelity Tr(rho @ GHZ) ', pformat(around(trace(GHZ_circuit_unitary @ GHZ_ket), 3)))
