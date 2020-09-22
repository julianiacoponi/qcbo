#!/usr/bin python3

import numpy as np
from qutip.qip.circuit import QubitCircuit, Gate
from qutip.qip.operations.gates import gate_sequence_product
from qutip import measurement

# NOTE: qip's circuit.py seems to do nothing meaningful with input or output states...
INPUT_STATES = ['0', '0', '0']
np.random.seed(28)

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

def create_fig_3_GHZ_circuit(phis=None):
    '''
    create the circuit in Fig. 3 of
    https://arxiv.org/pdf/1909.01229.pdf
    '''
    if phis is None:
        phis = np.random.uniform(0, 2 * np.pi, 6)

    qubit_circuit = QubitCircuit(3, input_states=INPUT_STATES)

    print('phis = ', phis)

    qubit_circuit.add_gate('RX', targets=0, arg_value=phis[0])
    qubit_circuit.add_gate('RX', targets=1, arg_value=phis[1])
    qubit_circuit.add_gate('RY', targets=2, arg_value=phis[2])

    qubit_circuit.add_gate('CNOT', controls=0, targets=1)
    qubit_circuit.add_gate('CNOT', controls=0, targets=2)

    qubit_circuit.add_gate('RX', targets=0, arg_value=phis[3])
    qubit_circuit.add_gate('RX', targets=1, arg_value=phis[4])
    qubit_circuit.add_gate('RY', targets=2, arg_value=phis[5])

    return qubit_circuit

def run_qubit_circuit(qubit_circuit):
    ''' run the given qubit circuit '''
    operator_matrix_list = qubit_circuit.propagators()
    return gate_sequence_product(operator_matrix_list)

if __name__ == '__main__':
    print('Running IBM Tutorial GHZ circuit ====================================')
    ibm_rho = run_qubit_circuit(create_IBM_GHZ_circuit())
    print(ibm_rho)

    print('Running Fig. 3 GHZ circuit ==========================================')
    fig3_rho = run_qubit_circuit(create_fig_3_GHZ_circuit())
    print(fig3_rho)
