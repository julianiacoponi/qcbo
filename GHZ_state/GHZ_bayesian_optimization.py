#!/usr/bin python3

import time
import numpy as np
from functools import partial
from itertools import count
from pprint import pformat

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
        # TODO: not sure about these bounds
        self.bounds = [(0, 1)] * self.input_dim
        self.min = [(0.) * self.input_dim]
        self.fmin = 0
        self.sd = standard_deviation

    def f_k(self, k, iteration_counter, X):
        # TODO: figure out how GPyOpt chooses the input here
        # something to do with Design_space class?
        # X = input phis
        X = reshape(X, self.input_dim)
        # run the GHZ circuit with the input phis to get the 8x8 triple qubit state matrix, rho
        rho = run_qubit_circuit(create_fig_3_GHZ_circuit(X))
        # evaluate the surrogate model
        fval = np.trace(rho * P_func(k))

        iteration_count = next(iteration_counter)
        # only print the first and then every % nth iteration
        if iteration_count == 1 or not iteration_count % 10:
            print(
                '\n' + '-' * 50 + f' iteration #{iteration_count}'
                f'\nGHZ f_{k}: X = phi_j column vector = {X}'

                f'\nGHZ f_{k}: phis (in degrees)'
                f' = {[np.round(phi * (180 / np.pi), 1) for phi in X[0]]}'

                f'\nGHZ f_{k}: fval = f(phis) = {fval}'
            )

        # TODO: do I want to add noise?
        n = X.shape[0]
        if self.sd == 0:
            noise = np.zeros(n).reshape(n, 1)
        else:
            noise = np.random.normal(0, self.sd, n).reshape(n, 1)
        return fval.reshape(n, 1)

if __name__ == '__main__':
    domain = [
        {'name': f'phi_{j}', 'type': 'continuous', 'domain': (0, 2 * np.pi)} for j in range(6)
    ]
    func = GHZ()
    surrogate_bopts = []

    for k in range(1, 8):
        bopt = GPyOpt.methods.BayesianOptimization(
            f=partial(func.f_k, k, count(1)),
            domain=domain,

            # aiming to maximize the first 4 surrogate models while minimising the last 3
            maximize=True if k <= 4 else False,

            # TODO: figure out if this is the right acquisition
            acquisiton_type='EI',

            # TODO: figure out if I want exact_feval (default is False)
            exact_feval=True,

            # NOTE: kernel choice by default seems to be GPy.kern.Matern52 as desired
            # this is the R_2(x) kernel of Eq.22 as applied to Eq.19 (of the paper)
        )

        surrogate_bopts.append(bopt)
        bopt.run_optimization(max_iter=100, eps=1e-6, max_time=300)
        print('\n\n' + '=' * 100)
        param_names = bopt.model.get_model_parameters_names()
        params = bopt.model.get_model_parameters()
        print(
            f'Optimal phis for model {k} = {bopt.x_opt}'
            f'\nOptimal f(phis) for model {k} = {bopt.fx_opt}'
            # this just shows final 6 phis?
            # f'\nphis for model {k} = {np.round(bopt.X, 2)}'
            f'\nmodel params names = {param_names}, len = {len(param_names)}'
            f'\nmodel params = {params}, len = {len(params[0])}'
        )
        print('=' * 100)

        bopt.plot_convergence(filename=f'/qcbo/k{k}_output_convergence.png')
        bopt.save_report(report_file=f'/qcbo/k{k}_output_report.txt')
        print('\n😴😴😴\nsleeping for 1s before next loop...\n😴😴😴\n')
        time.sleep(1)

    f_1, f_2, f_3, f_4, f_5, f_6, f_7 = tuple([bo.fx_opt for bo in surrogate_bopts]                                                                                                                                                                               )

    Fidelity = 1 / 4 * ((f_1 + f_2 + f_3 + f_4) - (f_5 + f_6 + f_7))
    print('Fidelity = ', Fidelity)
    print('Infidelity = ', 1 - Fidelity)
