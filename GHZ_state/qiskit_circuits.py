#!/usr/bin python3

from qiskit import QuantumCircuit, Aer, execute

def get_unitary_from_circuit(quantum_circuit, matrix_debug=False):
    ''' gets the unitary matrix for the given circuit '''
    backend = Aer.get_backend('unitary_simulator')
    return execute(quantum_circuit, backend).result().get_unitary(quantum_circuit)

class Bell:
    """ Class with the different circuits to produce the Bell state """

    def ideal(self, **kwargs):
        '''
        create the circuit which produce the Bell state
             ┌───┐
        q_0: ┤ H ├──■──
             └───┘┌─┴─┐
        q_1: ─────┤ X ├
                  └───┘
        '''
        qc = QuantumCircuit(2)
        qc.h(0)
        qc.cx(0, 1)

        if kwargs.get('debug'):
            print(f'QuantumCircuit:\n{qc.draw()}')

        return qc

    def rotations(self, angles, for_gypopt=True, **kwargs):
        '''
        Circuit to produce a Bell state with rotation gates
        `angles` should have a list of angles as its first element e.g. [[a1, a2, ..., a6]]
        (for compatability with GPyOpt)

        create the circuit in Fig. 3 of
        https://arxiv.org/pdf/1909.01229.pdf
        e.g. a perfect fidelity state is given by
             ┌─────────┐      ┌─────────┐
        q_0: ┤ RX(π/2) ├───■──┤ RX(π/2) ├
             ├─────────┴┐┌─┴─┐├─────────┤
        q_1: ┤ RX(3π/2) ├┤ X ├┤ RY(π/2) ├
             └──────────┘└───┘└─────────┘
        '''

        # this is to handle the input angles when using GPyOpt which are a 2d array of dimensions 1,6
        if for_gypopt:
            angle_1, angle_2, angle_3, angle_4 = angles[0]
        else:
            angle_1, angle_2, angle_3, angle_4 = angles

        qc = QuantumCircuit(2)

        qc.rx(angle_1, 0)
        qc.rx(angle_2, 1)
        qc.cx(0, 1)
        qc.rx(angle_3, 0)
        qc.ry(angle_4, 1)

        if kwargs.get('debug'):
            print(f'QuantumCircuit:\n{qc.draw()}')

        return qc


class GHZ:
    """ Class with the different circuits to produce the GHZ state """

    def ideal(self, **kwargs):
        '''
        create the circuit which produces a GHZ state
             ┌───┐
        q_0: ┤ H ├──■────■──
             └───┘┌─┴─┐  │
        q_1: ─────┤ X ├──┼──
                  └───┘┌─┴─┐
        q_2: ──────────┤ X ├
                       └───┘
        '''
        qc = QuantumCircuit(3)

        qc.h(0)
        qc.cx(0, 1)
        qc.cx(0, 2)

        if kwargs.get('debug'):
            print(f'QuantumCircuit:\n{qc.draw()}')

        return qc

    def ideal_2(self, **kwargs):
        '''
        create alternative equivalent ideal GHZ circuit
             ┌───┐┌───┐┌───┐
        q_0: ┤ H ├┤ X ├┤ H ├─────
             ├───┤└─┬─┘├───┤┌───┐
        q_1: ┤ H ├──┼──┤ X ├┤ H ├
             ├───┤  │  └─┬─┘├───┤
        q_2: ┤ H ├──■────■──┤ H ├
             └───┘          └───┘
        '''
        qc = QuantumCircuit(3)

        qc.h(0)
        qc.h(1)
        qc.h(2)
        qc.cx(2, 0)
        qc.cx(2, 1)
        qc.h(0)
        qc.h(1)
        qc.h(2)

        if kwargs.get('debug'):
            print(f'QuantumCircuit:\n{qc.draw()}')

        return qc

    def rotations(self, angles, for_gypopt=True, **kwargs):
        '''
        Circuit to produce a GHZ state with rotation gates
        `angles` should have a list of angles as its first element e.g. [[a1, a2, ..., a6]]
        (for compatability with GPyOpt)

        create the circuit in Fig. 3 of
        https://arxiv.org/pdf/1909.01229.pdf

        e.g. these angles which gives state (1/sqrt(2)) * (|000> - i*|111>)
             ┌─────────┐          ┌───────┐
        q_0: ┤ RX(π/2) ├──■────■──┤ RX(0) ├
             └┬───────┬┘┌─┴─┐  │  ├───────┤
        q_1: ─┤ RX(0) ├─┤ X ├──┼──┤ RX(0) ├
              ├───────┤ └───┘┌─┴─┐├───────┤
        q_2: ─┤ RX(0) ├──────┤ X ├┤ RY(0) ├
              └───────┘      └───┘└───────┘
        '''

        # this is to handle the input angles when using GPyOpt which are a 2d array of dimensions 1,6
        if for_gypopt:
            angle_1, angle_2, angle_3, angle_4, angle_5, angle_6 = angles[0]
        else:
            angle_1, angle_2, angle_3, angle_4, angle_5, angle_6 = angles

        qc = QuantumCircuit(3)

        qc.rx(angle_1, 0)
        qc.rx(angle_2, 1)
        qc.rx(angle_3, 2)
        qc.cx(0, 1)
        qc.cx(0, 2)
        qc.rx(angle_4, 0)
        qc.rx(angle_5, 1)
        qc.ry(angle_6, 2)

        if kwargs.get('debug'):
            print(f'QuantumCircuit:\n{qc.draw()}')

        return qc
