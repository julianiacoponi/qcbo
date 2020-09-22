#!/usr/bin python3
'''
Compute the Fidelity of a given 3-qubit state `rho` mapping to the GHZ state `Phi`
see Section III.B of: https://arxiv.org/pdf/1909.01229.pdf
GHZ state wikipedia: https://en.wikipedia.org/wiki/Greenberger%E2%80%93Horne%E2%80%93Zeilinger_state
'''
import numpy as np

# 2x2 identity matrix
I = np.eye(2)

# Pauli matrices
sx = np.array([[0, 1], [1, 0]])
sy = np.array([[0, -1j], [1j, 0]])
sz = np.array([[1, 0], [0, -1]])

zero = np.array([[1],[0]])  # aka 'spin up' or |0>
one = np.array([[0],[1]])  # aka 'spin down' or |1>


def triple_kronecker_product(matrix_1, matrix_2, matrix_3):
    ''' Does m1 (x) m2 (x) m3 '''
    return np.kron(np.kron(matrix_1, matrix_2), matrix_3)


triple_0 = triple_kronecker_product(zero, zero, zero)  # |000>
triple_1 = triple_kronecker_product(one, one, one)  # |111>


def S_func(k):
    '''
    components of the Fidelity equation for the GHZ state given in Eq. 26 of
    https://arxiv.org/pdf/1909.01229.pdf
    '''
    if k == 1:
        return triple_kronecker_product(sz, sz, sz)

    # NOTE: could use some itertools things here, but it's only 6 clauses so clearer to be explicit
    # 3 permutations of I (x) sz (x) sz
    if k == 2:
        return triple_kronecker_product(I, sz, sz)
    if k == 3:
        return triple_kronecker_product(sz, I, sz)
    if k == 4:
        return triple_kronecker_product(sz, sz, I)

    # 3 permutations of sx (x) sy (x) sy
    if k == 5:
        return triple_kronecker_product(sx, sy, sy)
    if k == 6:
        return triple_kronecker_product(sy, sx, sy)
    if k == 7:
        return triple_kronecker_product(sy, sy, sx)


def P_func(k):
    ''' measurement probability '''
    return (S_func(k) + I) / 2


def Fidelity(rho=None):
    ''' Fidelity for state rho resulting from GHZ gate sequence '''
    if rho is None:
        rho = triple_0
    positive_contributions = [np.trace(rho * S_func(k)) for k in range(1, 5)]
    negative_contributions = [np.trace(rho * S_func(k)) for k in range(5, 8)]
    return 1 / 8 * (1 + np.sum(positive_contributions) - np.sum(negative_contributions))


def Fidelity_in_terms_of_measurement_probablity(rho=None):
    if rho is None:
        rho = triple_0
    # TODO !!!


if __name__ == '__main__':
    max_k = 7
    S_list = [None  ]
    for k in range(1, max_k + 1):
        S_list.append(S_func(k))
        print(f'S_{k} = ', S_list[k])

    print('Fidelity for rho=|000> =', Fidelity())
    print('Fidelity for rho=|111> =', Fidelity(rho=triple_1))
