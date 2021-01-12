#!/usr/bin python3

from numpy import eye, array, kron, sqrt, sin, cos

# make i actually i :)
i = 1j

# 2x2, 4x4, 8x8 Identity matrices
I2, I4, I8 = eye(2), eye(4), eye(8)

#|0> or 'spin up'
ZERO = array([
    [1],
    [0],
])

#|1> or 'spin down'
ONE = array([
    [0],
    [1],
])

# Pauli matrices
Pauli_X = array([
    [0, 1],
    [1, 0],
])
Pauli_Y = array([
    [0, -i],
    [i, 0],
])
Pauli_Z = array([
    [1, 0],
    [0, -1],
])

HADAMARD = (1 / sqrt(2)) * array([[1, 1],[1, -1]])

CNOT = array([
    [1, 0, 0, 0],
    [0, 1, 0, 0],
    [0, 0, 0, 1],
    [0, 0, 1, 0],
])

def Rx(angle):
    ''' RX gate matrix '''
    a = angle / 2
    return array([
        [cos(a), -i*sin(a)],
        [-i*sin(a), cos(a)],
    ])


def Ry(angle):
    ''' RY gate matrix '''
    a = angle / 2
    return array([
        [cos(a), -sin(a)],
        [sin(a), cos(a)],
    ])

def triple_kron(matrix_1, matrix_2, matrix_3):
    ''' Does matrix_1 (x) matrix_2 (x) matrix_3 '''
    return kron(kron(matrix_1, matrix_2), matrix_3)


# CNOT_AB: if qubit A has value 1, flip qubit B
# if A = first qubit, this should be 4x4 Identity in the top left
# formula used: CNOT_01 = (|0><0| (x) I (x) I) + (|1><1| (x) Pauli_X (x) I)
# top left 4x4 corner of the 8x8 matrix
top_left = triple_kron(ZERO @ ZERO.conj().T, I2, I2)
# Pauli_X does the flipping
CNOT_01 = top_left + triple_kron(ONE @ ONE.conj().T, Pauli_X, I2)
CNOT_02 = top_left + triple_kron(ONE @ ONE.conj().T, I2, Pauli_X)


def qubit_state(qubit_string):
    '''
    Given a 2 or 3-string of 0s and 1s, provides the qubit column vector.
    This is all 0s except for a 1 at the binary value of the 2 or 3-string.
    e.g. 10 = 2 in binary, will have value 1 in the 3rd row (00 is 1 in the 1st row)
    e.g. 101 = 5 in binary, will have value 1 in the 6th row (000 is 1 in the 1st row)
    '''
    qbs = qubit_string

    if len(qbs) == 2:
        if qbs == '00':
            return kron(ZERO, ZERO)
        if qbs == '01':
            return kron(ZERO, ONE)
        if qbs == '10':
            return kron(ONE, ZERO)
        if qbs == '11':
            return kron(ONE, ONE)

    elif len(qbs) == 3:
        if qbs == '000':
            return triple_kron(ZERO, ZERO, ZERO)
        if qbs == '001':
            return triple_kron(ZERO, ZERO, ONE)
        if qbs == '010':
            return triple_kron(ZERO, ONE, ZERO)
        if qbs == '011':
            return triple_kron(ZERO, ONE, ONE)
        if qbs == '100':
            return triple_kron(ONE, ZERO, ZERO)
        if qbs == '101':
            return triple_kron(ONE, ZERO, ONE)
        if qbs == '110':
            return triple_kron(ONE, ONE, ZERO)
        if qbs == '111':
            return triple_kron(ONE, ONE, ONE)

    raise ValueError(f'`qubit_string` must be a length 2 or 3 string of 0s and 1s, not "{qbs}"')


BELL_STATE = (1 / sqrt(2)) * ((qubit_state('00') + qubit_state('11')))
GHZ_STATE = (1 / sqrt(2)) * (qubit_state('000') + qubit_state('111'))
