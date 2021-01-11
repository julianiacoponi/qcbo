#!/usr/bin python3
'''
Compute the Fidelity of either a given
- 2-qubit state mapping to the Bell state, or
- 3-qubit state mapping to the GHZ state

see Section III.B of: https://arxiv.org/pdf/1909.01229.pdf
GHZ state wikipedia: https://en.wikipedia.org/wiki/Greenberger%E2%80%93Horne%E2%80%93Zeilinger_state
'''
import time
import itertools
from math import isclose
from pprint import pformat
from argparse import ArgumentParser
from numpy import eye, array, kron, sqrt, sin, cos, pi, around, trace, random, linalg
from qutip.qip.circuit import QubitCircuit, Gate
from qutip.qip.operations.gates import gate_sequence_product


parser = ArgumentParser('GHZ_fidelity_with_numpy')
parser.add_argument(
    '--input_qubits',
    help=(
        '2 or 3-string of 0s and 1s for the 3 input qubits.'
        ' If 2 qubits, the Bell circuit is used instead'
    ),
    default=None,
)
parser.add_argument(
    '--bell',
    dest='bell',
    help='Calculate fidelity of the Bell circuit state (overrides --ghz)',
    action='store_true',
    default=False,
)
parser.add_argument(
    '--ghz',
    dest='ghz',
    help='Calculate fidelity of the GHZ circuit state (overridden by --bell)',
    action='store_true',
    default=False,
)
parser.add_argument(
    '--angles',
    nargs='*',
    help=(
        '4 or 6 angles (default units of π) for the Bell or GHZ circuit rotation gates'
        '. Random angles chosen'
    ),
)
parser.add_argument(
    '--units',
    dest='units',
    help='units for the input angles',
    default='pi',
)
parser.add_argument(
    '--permute',
    dest='permute',
    action='store_true',
    default=False,
    help=(
        'permute (with repititions) and iterate over the given `--angles`'
        '. If `--bell`, will permute a choice of 4. If `--ghz`, will permute a choice of 6'
        '. This means you can run the script with e.g. `--bell --angles 1 2 3 4 5 6 --permute`'
        ', and all 6^4 permutations of 4 choices (for bell) of 6 angles (provided)'
        ' will be iterated over.'
    ),
)
parser.add_argument(
    '--with_qutip',
    dest='with_qutip',
    help='use unitary generated with QubitCircuit from qutip',
    action='store_true',
    default=False,
)
parser.add_argument(
    '--debug',
    dest='debug',
    help='print extra debug',
    action='store_true',
    default=False,
)
parser.add_argument(
    '--matrix_debug',
    dest='matrix_debug',
    help='print extra matrix debug',
    action='store_true',
    default=False,
)
parser.add_argument(
    '--random_seed',
    help='argument for numpy.random.seed, to initiate choices for --random_angles',
    type=int,
    default=28,
)

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

def Bell_circuit_unitary(angles, with_qutip=False, matrix_debug=False):
    '''
    Bell circuit unitary
    `angles` has a list of angles as its first element e.g. [[a1, a2, a3, a4]]
    (for compatability with GPyOpt)

    if with_qutip:
    create the Bell state equivalent of the GHZ circuit in Fig. 3 of
    https://arxiv.org/pdf/1909.01229.pdf
    i.e.
    ---[Rx_1]---*---[Rx_3]---
                |
    ---[Rx_2]--(+)--[Ry_4]---
    '''
    assert_angles_in_radians(angles)
    angle_1, angle_2, angle_3, angle_4 = angles[0]

    if not with_qutip:
        Rx_1 = kron(Rx(angle_1), I2)
        Rx_2 = kron(I2, Rx(angle_2))
        Rx_3 = kron(Rx(angle_3), I2)
        Ry_4 = kron(I2, Ry(angle_4))

        if matrix_debug:
            print('\n')
            print(f'Rx_1 =\n{around(Rx_1, 3)}')
            print(f'Rx_2 =\n{around(Rx_2, 3)}')
            print(f'CNOT=\n{around(CNOT, 3)}')
            print(f'Rx_3 =\n{around(Rx_3, 3)}')
            print(f'Ry_4 =\n{around(Ry_4, 3)}')

        # apply the gates in sequence from the left
        return Ry_4 @ Rx_3 @ CNOT @ Rx_2 @ Rx_1

    if debug:
        print('-' * 100)
        print('Getting the Bell circuit unitary with qutip')

    qubit_circuit = QubitCircuit(2)

    angle_1, angle_2, angle_3, angle_4 = angles[0]

    qubit_circuit.add_gate('RX', targets=0, arg_value=angle_1)
    qubit_circuit.add_gate('RX', targets=1, arg_value=angle_2)

    qubit_circuit.add_gate('CNOT', controls=0, targets=1)

    qubit_circuit.add_gate('RX', targets=0, arg_value=angle_3)
    qubit_circuit.add_gate('RY', targets=1, arg_value=angle_4)

    operator_matrix_list = qubit_circuit.propagators()
    if matrix_debug:
        print('operator matrix list = \n', pformat(operator_matrix_list))

    return gate_sequence_product(operator_matrix_list)


def GHZ_circuit_unitary(angles, with_qutip=False, matrix_debug=False):
    '''
    GHZ circuit unitary
    `angles` has a list of angles as its first element e.g. [[a1, a2, ..., a6]]
    (for compatability with GPyOpt)

    if with_qutip:
    create the circuit in Fig. 3 of
    https://arxiv.org/pdf/1909.01229.pdf
    ---[Rx_1]---*---*---[Rx_4]---
                |   |
    ---[Rx_2]--(+)--|---[Rx_5]---
                    |
    ---[Rx_3]------(+)--[Ry_6]---

    '''
    assert_angles_in_radians(angles)
    angle_1, angle_2, angle_3, angle_4, angle_5, angle_6 = angles[0]

    if not with_qutip:
        Rx_1 = triple_kron(Rx(angle_1), I2, I2)
        Rx_2 = triple_kron(I2, Rx(angle_2), I2)
        Rx_3 = triple_kron(I2, I2, Rx(angle_3))

        Rx_4 = triple_kron(Rx(angle_4), I2, I2)
        Rx_5 = triple_kron(I2, Rx(angle_5), I2)
        Ry_6 = triple_kron(I2, I2, Ry(angle_6))

        if matrix_debug:
            print('\n')
            print(f'Rx_1 =\n{around(Rx_1, 3)}')
            print(f'Rx_2 =\n{around(Rx_2, 3)}')
            print(f'Rx_3 =\n{around(Rx_3, 3)}')
            print(f'CNOT_01 =\n{around(CNOT_01, 3)}')
            print(f'CNOT_2 =\n{around(CNOT_02, 3)}')
            print(f'Rx_4 =\n{around(Rx_4, 3)}')
            print(f'Rx_5 =\n{around(Rx_5, 3)}')
            print(f'Ry_6 =\n{around(Ry_6, 3)}')

        # apply the gates in sequence from the left
        return Ry_6 @ Rx_5 @ Rx_4 @ CNOT_02 @ CNOT_01 @ Rx_3 @ Rx_2 @ Rx_1

    if debug:
        print('-' * 100)
        print('Getting the GHZ circuit unitary with qutip')

    qubit_circuit = QubitCircuit(3)

    angle_1, angle_2, angle_3, angle_4, angle_5, angle_6 = angles[0]

    qubit_circuit.add_gate('RX', targets=0, arg_value=angle_1)
    qubit_circuit.add_gate('RX', targets=1, arg_value=angle_2)
    qubit_circuit.add_gate('RX', targets=2, arg_value=angle_3)

    qubit_circuit.add_gate('CNOT', controls=0, targets=1)
    qubit_circuit.add_gate('CNOT', controls=0, targets=2)

    qubit_circuit.add_gate('RX', targets=0, arg_value=angle_4)
    qubit_circuit.add_gate('RX', targets=1, arg_value=angle_5)
    qubit_circuit.add_gate('RY', targets=2, arg_value=angle_6)

    operator_matrix_list = qubit_circuit.propagators()
    if matrix_debug:
        print('operator matrix list = \n', pformat(operator_matrix_list))

    return gate_sequence_product(operator_matrix_list)


def assert_angles_in_radians(angles):
    ''' Sanity check that angles is in radians '''
    for angle in angles[0]:
        assert 0 <= angle <= 2 * pi, f'{around(angle, 3)} is negative or more than 2π, are you sure it is in radians?'


def density_matrix_from_unitary(unitary, input_qubits, with_qutip=False):
    '''
    Calculate the density matrix for a given circuit unitary and input qubit state
    Density matrix rho = U|psi><psi|U`
    where U is the circuit unitary, U` = conjugate transpose of U
    |psi> is the input (triple) qubit state to the circuit producing the unitary
    '''
    input_qubit_state = qubit_state(qubit_string=input_qubits)
    if with_qutip:
        # convert the Qobj's into arrays first
        return array(unitary) @ input_qubit_state @ input_qubit_state.conj().T @ array(unitary.conj().trans())

    return unitary @ input_qubit_state @ input_qubit_state.conj().T @ unitary.conj().T


def S_operator(k):
    '''
    components of the Fidelity equation for the GHZ state given in Eq. 26 of
    https://arxiv.org/pdf/1909.01229.pdf
    '''
    if k == 1:
        return triple_kron(Pauli_X, Pauli_X, Pauli_X)

    # 3 permutations of I (x) Pauli_Z (x) Pauli_Z
    if k == 2:
        return triple_kron(I2, Pauli_Z, Pauli_Z)
    if k == 3:
        return triple_kron(Pauli_Z, I2, Pauli_Z)
    if k == 4:
        return triple_kron(Pauli_Z, Pauli_Z, I2)

    # 3 permutations of Pauli_X (x) Pauli_Y (x) Pauli_Y
    if k == 5:
        return triple_kron(Pauli_X, Pauli_Y, Pauli_Y)
    if k == 6:
        return triple_kron(Pauli_Y, Pauli_X, Pauli_Y)
    if k == 7:
        return triple_kron(Pauli_Y, Pauli_Y, Pauli_X)


def P_operator(k):
    ''' measurement probability of positive value S_operator'''
    # adding the identity then halving eliminates all the -1s on the diagonal, leaving only +1s
    return (S_operator(k) + I8) / 2


def Fidelity_Sk(rho):
    ''' Fidelity for state rho resulting from GHZ gate sequence '''
    positive_contributions = [trace(rho @ S_operator(k)) for k in range(1, 5)]
    negative_contributions = [trace(rho @ S_operator(k)) for k in range(5, 8)]
    return 1 / 8 * (1 + sum(positive_contributions) - sum(negative_contributions))


def Fidelity_Pk(rho, debug=False):
    ''' Fidelity for state rho resulting from GHZ gate sequence '''
    f1, f2, f3, f4 = [trace(rho @ P_operator(k)) for k in range(1, 5)]
    f5, f6, f7 = [trace(rho @ P_operator(k)) for k in range(5, 8)]
    probablities = [f1, f2, f3, f4, f5, f6, f7]

    if debug:
        for k, prob in enumerate(probablities, 1):
            print(f'f{k}={around(prob, 3)}')


    return 1 / 4 * ((f1 + f2 + f3 + f4) - (f5 + f6 + f7))

def convert_angles(angles, unit, round_dp=3, for_gpyopt=True):
    ''' Converts angles in one unit to the other two, of π,rad,º, returned in that order'''
    if unit == 'pi':
        angles_pi = angles
        angles_radians = [angle * pi for angle in angles]
        angles_degrees = [angle * 180 for angle in angles]

    elif unit == 'radians':
        angles_pi = [angle / pi for angle in angles]
        angles_radians = angles
        angles_degrees = [angle * (180 / pi) for angle in angles]

    elif unit == 'degrees':
        angles_pi = [angle / 180 for angle in angles]
        angles_radians = [angle * (pi / 180) for angle in angles]
        angles_degrees = angles

    else:
        raise AttributeError(f'`unit` must be "pi", "radians" or "degrees", not "{unit}"')

    # GPyOpt requires [[angle_1, angle_2, ...]]
    if for_gpyopt:
        return [angles_pi], [angles_radians], [angles_degrees]

    return angles_pi, angles_radians, angles_degrees

def calculate_state_fidelity(
    bell,
    ghz,
    input_qubits,
    angles,
    units,
    with_qutip=False,
    debug=False,
    matrix_debug=False,
):
    ''' the main script function '''
    if bell and ghz:
        print(
            'Warning: Both Bell and GHZ requested'
            ', but Bell takes precedent'
            ', so not calculating GHZ'
        )
    if bell:
        assert len(input_qubits) == 2, 'Bell circuit needs 2-qubit input state'
        assert len(angles) == 4, 'Bell circuit needs 4 angles'
        ghz = False
    elif ghz:
        assert len(input_qubits) == 3, 'GHZ circuit needs 3-qubit input state'
        assert len(angles) == 6, 'GHZ circuit needs 6 angles'

    angles_pi, angles_radians, angles_degrees = convert_angles(angles, units)
    if debug:
        print('-' * 100)
        print(f'Input angles/π) = {around(angles_pi, 3)}')
        print(f'Input angles/radians) = {around(angles_radians, 3)}')
        print(f'Input angles/º) = {around(angles_degrees, 1)}')

    if bell:
        state_str = 'Bell'
        desired_state = BELL_STATE
        desired_state_str = '|00> + |11>'
        # ideal Bell circuit
        # |0> ---[H]---*---
        #              |
        # |0> --------(+)--
        ideal_unitary = CNOT @ kron(HADAMARD, I2)
        Unitary = Bell_circuit_unitary(
            angles_radians,
            with_qutip=with_qutip,
            matrix_debug=matrix_debug,
        )
    elif ghz:
        state_str = 'GHZ'
        desired_state = GHZ_STATE
        desired_state_str = '|000> + |111>'
        # ideal GHZ circuit
        # |0> ---[H]---*---*---
        #              |   |
        # |0> --------(+)--|---
        #                  |
        # |0> ------------(+)--
        ideal_unitary = CNOT_02 @ CNOT_01 @ triple_kron(HADAMARD, I2, I2)
        Unitary = GHZ_circuit_unitary(
            angles_radians,
            with_qutip=with_qutip,
            matrix_debug=matrix_debug,
        )

    input_state = qubit_state(input_qubits)
    if matrix_debug:
        print(f'\nIdeal {state_str} circuit unitary * sqrt(2) =\n', around(sqrt(2) * ideal_unitary, 3))
        print(f'\n{state_str} circuit unitary =\n{around(Unitary, 3)}')

        print(f'\nInput qubit state = |{input_qubits}> =\n', input_state)
        print(
            f'\nDesired |{state_str}> * sqrt(2)'
            f' = {desired_state_str} =\n',
            around(sqrt(2) * desired_state, 3),
        )

    output_state = Unitary @ input_state
    rho = density_matrix_from_unitary(Unitary, input_qubits, with_qutip=with_qutip)
    if debug:
        print(f'\nOutput qubit state = Unitary.|{input_qubits}> =\n', around(output_state, 3))
        print(f'\n{state_str} density matrix (rho) =\n{around(rho, 3)}')

    # Sk and Pk only apply to the GHZ state
    if ghz:
        F_Sk, F_Pk = Fidelity_Sk(rho), Fidelity_Pk(rho)
        if matrix_debug:
            print('\n')
            for k in range(1, 8):
                print(f'S_{k} = \n{S_operator(k)}')
                print(f'S_{k} + I =\n{S_operator(k) + I8}')
                print(f'P_{k} = (S_{k} + I)/2\n{P_operator(k)}')

            # off by ~2e-17
            print(f'before rounding {F_Sk=}')
            print(f'before rounding {F_Pk=}')

        assert F_Sk.imag == 0
        assert F_Pk.imag == 0
        if debug:
            print('Fidelity with Sk =', around(F_Sk, 3))
            print('Fidelity with Pk =', around(F_Pk, 3))

    Fidelity_matrix = desired_state.conj().T @ rho @ desired_state

    # Sanity checks, there should just be 1 real Fidelity value calculated!
    assert Fidelity_matrix.size == 1
    assert len(Fidelity_matrix[0]) == 1
    assert Fidelity_matrix[0][0].imag == 0

    Fidelity = Fidelity_matrix[0][0].real
    if debug:
        print(f'Fidelity <{state_str}|rho|{state_str}> =', around(Fidelity, 3))

    if matrix_debug:
        print('Is rho really a density matrix? Tr(rho) == 1 should be True, is', around(trace(rho), 2) == 1)
        trace_rho_squared = around(trace(rho @ rho), 2)
        is_pure = trace_rho_squared == 1
        is_mixed = trace_rho_squared < 1
        print(
            f'Does rho represent a pure or mixed state?'
            f' Tr(rho^2) == 1 is pure, < 1 mixed'
            f'; is {trace_rho_squared}, so pure={is_pure}, mixed={is_mixed}'
        )
    return Fidelity, output_state


if __name__ == '__main__':
    t0 = time.perf_counter()
    args = parser.parse_args()
    print('-' * 100)
    print('Command line args are', args)
    print('-' * 100)

    # either or (bell takes precedent if both are set)
    bell = args.bell
    ghz = args.ghz

    input_qubits = args.input_qubits  # default None
    angles = args.angles  # should be length 4 for bell, 6 for ghz
    units = args.units  # default 'pi'
    permute = args.permute  # default False

    with_qutip = args.with_qutip  # default False

    random_seed = args.random_seed  # default 28
    if random_seed != 0:
        random.seed(random_seed)

    debug = args.debug
    matrix_debug = args.matrix_debug

    if not input_qubits:
        if bell:
            input_qubits = '00'
        elif ghz:
            input_qubits = '000'

    num_angles = 4 if bell else 6
    if not angles:
        limit = {'pi': 2, 'radians': 2 * pi, 'degrees': 360}
        if bell:
            angles = list(random.uniform(0, limit[units], num_angles))
        elif ghz:
            angles = list(random.uniform(0, limit[units], num_angles))
    else:
        angles = [float(angle) for angle in angles]

    if permute:
        tic = time.perf_counter()
        print(f'Took {tic - t0:0.4f} seconds to start permutation loop')
        print(
            f'Starting iteration through'
            f' {len(angles)}**{num_angles}={len(angles)**num_angles} sets of angles...'
        )
        all_fidelities = {}
        all_output_states = {}

        # this means form 'ab' you get 'aa' 'ab' 'ba' 'bb', instead of just permutations 'ab' 'ba'
        angle_perms_with_repititions = itertools.product(angles, repeat=num_angles)
        count = 0
        for perm in angle_perms_with_repititions:
            count += 1
            if not count % 1000:
                print(f'On permutation #{count}')
            if debug:
                print('\n')
                print('_' * 100)
                print(f'PERMUTATION #{count}')

            angles_choice = list(perm)
            fidelity, output_state = calculate_state_fidelity(
                bell,
                ghz,
                input_qubits,
                angles_choice,
                units,
                with_qutip=with_qutip,
                debug=debug,
                matrix_debug=matrix_debug,
            )
            all_fidelities[str(angles_choice)] = round(fidelity, 3)
            all_output_states[str(angles_choice)] = around(output_state, 3)

        toc = time.perf_counter()
        print(
            f'{len(angles)}**{num_angles} = {count} permutations of {angles}/{units}'
            f' done in {toc - tic:0.4f} seconds'
        )

        # only print this if there's not thaaaat much output
        if count < 1000:
            print(f'\nAll fidelities for each set of angles:\n{pformat(all_fidelities)}')
            ordered_by_fidelity = [
                (angle_choice, all_fidelities[angle_choice])
                for angle_choice in sorted(
                    all_fidelities,
                    key=lambda ac: all_fidelities[ac],
                    reverse=True,
                )
            ]
            print(f'\nAll fidelities ordered:\n{pformat(ordered_by_fidelity)}')

        # NOTE: Found 4 perfect matches for Bell state, with angles 0, π/2, π, 3π/2
        perfect_fidelities = dict(filter(
                lambda key_value_pair: key_value_pair[1] == 1,
                all_fidelities.items(),
            ))
        print(
            f'\nAll {len(perfect_fidelities)} sets of angles with Fidelity == 1:'
            f'\n{pformat(perfect_fidelities)}'
        )

        half_perfect_fidelities = dict(filter(
                lambda key_value_pair: 0.5 <= key_value_pair[1] < 1,
                all_fidelities.items(),
            ))
        print(
            f'\nAll {len(half_perfect_fidelities)} sets of angles with 0.5 <= Fidelity < 1:'
            f'\n{pformat(half_perfect_fidelities)}'
        )

        if ghz:
            # want to find the angles that give just |000> and |111> (but with phase difference)
            # filter by ensuring all 6 other states have 0 in the state vector
            close_states = dict(filter(
                lambda key_value_pair: all(
                    [key_value_pair[1][n] == 0 for n in range(2, 7)]
                    + [key_value_pair[1][n] != 0 for n in (0, 7)]

                    # TODO: only pick states which have real or imaginary parts == 1/sqrt(2)
                    # + [
                    #     any([isclose(key_value_pair[1][n].real**2, 0.5, abs_tol=1e-4) for n in (0, 7)]),
                    #     any([isclose(key_value_pair[1][n].imag**2, 0.5, abs_tol=1e-4) for n in (0, 7)]),
                    # ]

                ),
                all_output_states.items(),
            ))

            print(
                f'\nAll {len(close_states)} sets of angles that give a mixed state of only |000> and |111>:'
                f'\n{pformat(close_states)}'
            )
            print(
                f'\nAll {len(close_states)} sets of angles'
                f' that give a mixed state of only |000> and |111> (angles only):'
            )
            for angles_choice, state in close_states.items():
                print(angles_choice)
            # NOTE: Found 32 close matches for GHZ state, with angles 0, π/2, π, 3π/2
            # annoyingly, all have a phase difference between the two states!

    else:
        fidelity, output_state = calculate_state_fidelity(
            bell,
            ghz,
            input_qubits,
            angles,
            units,
            with_qutip=with_qutip,
            debug=debug,
            matrix_debug=matrix_debug,
        )
