#!/usr/bin python3

import numpy as np
from functools import partial

import GPyOpt
from GPyOpt.models import BOModel
from GPyOpt.acquisitions.base import AcquisitionBase
from GPyOpt.core.task.cost import constant_cost_withGradients

from GHZ_circuit_with_qutip import create_fig_3_GHZ_circuit, run_qubit_circuit
from GHZ_fidelity_with_numpy import P_func

# the Fidelity in the GHZ state example is
# F(theta) = 1/4 (sum1,4 f_k(theta) - sum5,7 f_k(theta))
# the 'surrogate models' in GHZ state example are
# f_k(theta) = trace(rho(theta_j=1,6) * P_func(k)) ===> k=1 to 7, these are my 7 surrogate models

def reshape(x, input_dim):
    ''' Reshapes x into a matrix with input_dim columns '''
    x = np.array(x)
    if x.size == input_dim:
        x = x.reshape((1, input_dim))
    return x

class GHZ:
    ''' GHZ state fidelity '''
    def __init__(self, input_dim=8, standard_deviation=0):
        self.input_dim = input_dim
        self.bounds = [(0, 1)] * self.input_dim
        self.min = [(0.) * self.input_dim]
        self.fmin = 0
        self.sd = standard_deviation

    def f_k(self, k, X):
        # X = input phis
        X = reshape(X, self.input_dim)
        print('\n' + '-' * 50)
        print(f'GHZ f_{k}: X = {X}')
        n = X.shape[0]
        fval = np.trace(run_qubit_circuit(create_fig_3_GHZ_circuit(X)) * P_func(k))
        print(f'GHZ f_{k}: fval = {fval}')

        # TODO: do I want to add noise?
        if self.sd == 0:
            noise = np.zeros(n).reshape(n, 1)
        else:
            noise = np.random.normal(0, self.sd, n).reshape(n, 1)
        return fval.reshape(n, 1)

if __name__ == '__main__':

    # TODO: figure out if this domain is correct?
    domain = [
        {'name': f'phi_{j}', 'type': 'continuous', 'domain': (0, 2 * np.pi)} for j in range(6)
    ]

    func = GHZ()
    surrogate_bopts = []
    for k in range(1, 8):
        bopt = GPyOpt.methods.BayesianOptimization(
            f=partial(func.f_k, k),
            domain=domain,
            maximize=True,
            # TODO: figure out of these are the right ones
            acquisiton_type='EI',
            exact_feval=True,
        )
        surrogate_bopts.append(bopt)

    for k, bopt in enumerate(surrogate_bopts, 1):
        bopt.run_optimization(max_iter=20, eps=1e-6, max_time=3)
        print('\n\n' + '=' * 50)
        print(f'Optimal phi for model {k} = {bopt.x_opt}')
        print(f'phis for model {k} = {np.round(bopt.X, 2)}')
        bopt.plot_convergence(filename=f'/qcbo/k{k}_output_convergence.png')
        bopt.save_report(report_file=f'/qcbo/k{k}_output_report.txt')
