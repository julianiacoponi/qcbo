#!/usr/bin python3

from pprint import pformat
from numpy import random, pi, around
from qutip.qip.circuit import QubitCircuit, Gate
from qutip.qip.operations.gates import gate_sequence_product

# NOTE: qip's circuit.py seems to do nothing meaningful with input or output states...
INPUT_STATES = ['0', '0', '0']
random.seed(28)

# TODO: unsure how (or whether) to `add_measurement` to these circuits... probably?


def create_IBM_GHZ_circuit():
    '''
    create the circuit found at
    https://quantum-computing.ibm.com/docs/guide/mult-entang/ghz-states
    '''
    qubit_circuit = QubitCircuit(3, input_states=INPUT_STATES)

    qubit_circuit.add_gate('SNOT', targets=0)
    qubit_circuit.add_gate('SNOT', targets=1)
    qubit_circuit.add_gate('X', targets=2)

    qubit_circuit.add_gate('CNOT', controls=2, targets=0)
    qubit_circuit.add_gate('CNOT', controls=2, targets=1)

    qubit_circuit.add_gate('SNOT', targets=0)
    qubit_circuit.add_gate('SNOT', targets=1)
    qubit_circuit.add_gate('SNOT', targets=2)

    return qubit_circuit


def get_random_angles(n):
    ''' create n random angles '''
    return random.uniform(0, 2 * pi, n)


def create_fig_3_GHZ_circuit(phis=None):
    '''
    create the circuit in Fig. 3 of
    https://arxiv.org/pdf/1909.01229.pdf
    '''
    if phis is None:
        phis = get_random_angles(6)

    # this is to handle the input phis when using GPyOpt which are a 2d array of dimensions 1,6
    if getattr(phis, 'shape'):
        if phis.shape[0] == 1:
            phis = phis[0]

    qubit_circuit = QubitCircuit(3, input_states=INPUT_STATES)

    # phis = [pi / 3] * 6
    # phis = [pi, pi/2, pi/3, pi/4, pi/5, pi/6]
    phis_deg = [69.1, 220.1, 131.7, 160.3, 188.8, 80.9]
    phis = [phi * (pi/180) for phi in phis_deg]
    # don't print when doing large numbers of iterations
    print('phis (in degrees) = ', [round(phi * (180 / pi), 1) for phi in phis])

    qubit_circuit.add_gate('RX', targets=0, arg_value=phis[0])
    qubit_circuit.add_gate('RX', targets=1, arg_value=phis[1])
    qubit_circuit.add_gate('RX', targets=2, arg_value=phis[2])

    qubit_circuit.add_gate('CNOT', controls=0, targets=1)
    qubit_circuit.add_gate('CNOT', controls=0, targets=2)

    qubit_circuit.add_gate('RX', targets=0, arg_value=phis[3])
    qubit_circuit.add_gate('RX', targets=1, arg_value=phis[4])
    qubit_circuit.add_gate('RY', targets=2, arg_value=phis[5])

    return qubit_circuit


def run_qubit_circuit(qubit_circuit):
    ''' run the given qubit circuit '''
    operator_matrix_list = qubit_circuit.propagators()
    # print('operator matrix list', pformat(operator_matrix_list))
    return gate_sequence_product(operator_matrix_list)


if __name__ == '__main__':
    # print('Running IBM Tutorial GHZ circuit ====================================')
    # ibm_rho = run_qubit_circuit(create_IBM_GHZ_circuit())
    # print(ibm_rho)

    print('Running Fig. 3 GHZ circuit ==========================================')
    fig3_rho = run_qubit_circuit(create_fig_3_GHZ_circuit())
    # print(fig3_rho)
    print('|000> column:\n', around(fig3_rho[:,0], 3))
    print('|001> column:\n', around(fig3_rho[:,1], 3))
    print('|010> column:\n', around(fig3_rho[:,2], 3))
    print('|011> column:\n', around(fig3_rho[:,3], 3))
    print('|100> column:\n', around(fig3_rho[:,-4], 3))
    print('|101> column:\n', around(fig3_rho[:,-3], 3))
    print('|110> column:\n', around(fig3_rho[:,-2], 3))
    print('|111> column:\n', around(fig3_rho[:,-1], 3))
