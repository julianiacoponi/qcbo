#!/usr/bin python3

from pprint import pformat
from qutip.qip.circuit import QubitCircuit
from qutip.qip.operations.gates import gate_sequence_product


def get_unitary_from_circuit(qubit_circuit, **kwargs):
    ''' gets the unitary matrix for the given circuit '''
    operator_matrix_list = qubit_circuit.propagators()
    if kwargs.get('matrix_debug'):
        print('operator matrix list =\n', pformat(operator_matrix_list))
    return gate_sequence_product(operator_matrix_list)


class Bell:
    """ Class with the different circuits to produce the Bell state """

    @property
    def ideal(self):
        '''
        create the circuit which produce the Bell state
        |0> ---[H]---*---
                     |
        |0> --------(+)--
        '''
        qc = QubitCircuit(2)
        qc.add_gate('SNOT', targets=0)
        qc.add_gate('CNOT', controls=2, targets=0)
        qc.add_gate('CNOT', controls=2, targets=1)
        return qc

    def rotations(self, angles, for_gypopt=True):
        '''
        Circuit to produce a Bell state with rotation gates
        `angles` should have a list of angles as its first element e.g. [[a1, a2, ..., a6]]
        (for compatability with GPyOpt)

        create the circuit in Fig. 3 of
        https://arxiv.org/pdf/1909.01229.pdf
        ---[Rx_1]---*-----[Rx_3]---
                    |
        ---[Rx_2]--(+)----[Rx_4]---
        '''

        # this is to handle the input angles when using GPyOpt which are a 2d array of dimensions 1,6
        if for_gypopt:
            angle_1, angle_2, angle_3, angle_4 = angles[0]
        else:
            angle_1, angle_2, angle_3, angle_4 = angles

        qc = QubitCircuit(2)

        qc.add_gate('RX', targets=0, arg_value=angle_1)
        qc.add_gate('RX', targets=1, arg_value=angle_2)

        qc.add_gate('CNOT', controls=0, targets=1)

        qc.add_gate('RX', targets=1, arg_value=angle_3)
        qc.add_gate('RY', targets=2, arg_value=angle_4)

        return qc


class GHZ:
    """ Class with the different circuits to produce the GHZ state """

    @property
    def ideal(self):
        '''
        create the circuit which produces a GHZ state
        |0> --[H]---*---*---
                    |   |
        |0> -------(+)--|---
                        |
        |0> -----------(+)--
        '''
        qc = QubitCircuit(3)
        qc.add_gate('SNOT', targets=0)
        qc.add_gate('CNOT', controls=0, targets=1)
        qc.add_gate('CNOT', controls=0, targets=2)
        return qc

    @property
    def ideal_2(self):
        '''
        create alternative equivalent ideal GHZ circuit
        Note: SNOT == Hadamard
        |0> --[H]--(+)------[H]--
                    |
        |0> --[H]------(+)--[H]--
                    |   |
        |0> --[X]---*---*---[H]--
        '''
        qc = QubitCircuit(3)

        qc.add_gate('SNOT', targets=0)
        qc.add_gate('SNOT', targets=1)
        qc.add_gate('X', targets=2)

        qc.add_gate('CNOT', controls=2, targets=0)
        qc.add_gate('CNOT', controls=2, targets=1)

        qc.add_gate('SNOT', targets=0)
        qc.add_gate('SNOT', targets=1)
        qc.add_gate('SNOT', targets=2)

        return qc

    def rotations(self, angles, for_gypopt=True):
        '''
        Circuit to produce a GHZ state with rotation gates
        `angles` should have a list of angles as its first element e.g. [[a1, a2, ..., a6]]
        (for compatability with GPyOpt)

        create the circuit in Fig. 3 of
        https://arxiv.org/pdf/1909.01229.pdf
        ---[Rx_1]---*---*---[Rx_4]---
                    |   |
        ---[Rx_2]--(+)--|---[Rx_5]---
                        |
        ---[Rx_3]------(+)--[Ry_6]---
        '''

        # this is to handle the input angles when using GPyOpt which are a 2d array of dimensions 1,6
        if for_gypopt:
            angle_1, angle_2, angle_3, angle_4, angle_5, angle_6 = angles[0]
        else:
            angle_1, angle_2, angle_3, angle_4, angle_5, angle_6 = angles

        qc = QubitCircuit(3)

        qc.add_gate('RX', targets=0, arg_value=angle_1)
        qc.add_gate('RX', targets=1, arg_value=angle_2)
        qc.add_gate('RX', targets=2, arg_value=angle_3)

        qc.add_gate('CNOT', controls=0, targets=1)
        qc.add_gate('CNOT', controls=0, targets=2)

        qc.add_gate('RX', targets=0, arg_value=angle_4)
        qc.add_gate('RX', targets=1, arg_value=angle_5)
        qc.add_gate('RY', targets=2, arg_value=angle_6)

        return qc
