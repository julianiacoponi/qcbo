#!/usr/bin python3

from argparse import ArgumentParser
import os
import time
import numpy as np
from numpy import array, kron, eye, sqrt, cos, sin, trace, zeros, pi, random, around
from functools import partial
from itertools import count
from pprint import pformat

import GPyOpt

from GHZ_fidelity_with_numpy import (
    GHZ_circuit_unitary,
    density_matrix_from_unitary,
    S_operator,
    P_operator,
    Fidelity_Sk,
    Fidelity_Pk,
    convert_angles,
)

parser = ArgumentParser('GHZ_bayesian_optimization')
parser.add_argument(
    '--input_qubits',
    help='3-string of 0s and 1s for the 3 input qubits',
    default='000',
)
# TODO: Add BellStateFidelity class?
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
    default=True,
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
# TODO: remove this
parser.add_argument(
    '--separate',
    help='optimize the 7 surrogate models separately',
    action='store_true',
    default=False,
)
parser.add_argument(
    '--with_qutip',
    dest='with_qutip',
    help='use qutip QubitCircuit unitary from GHZ_circuit_with_qutip.py',
    action='store_true',
    default=False,
)
parser.add_argument(
    '--random_seed',
    help='argument for numpy.random.seed, to initiate choices for --random_angles',
    type=int,
    default=28,
)
parser.add_argument(
    '--print_every_nth_iteration',
    help='only print results for every nth iteration',
    type=int,
    default=5,
)
parser.add_argument(
    '--use_Sk',
    dest='use_Sk',
    help='use Sk for Fidelity calculation instead of Pk',
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
    '--max_iter',
    help='number of acquisitions the optimization goes through',
    type=int,
    default=5,
)

class GhzStateFidelity:
    '''
    6-dimensional bayesian optimization experiment class for calculating Eq. 28
    Modelled off experimentsNd.py
    '''

    def __init__(
        self,
        input_dim=6,
        standard_deviation=0,
        with_qutip=False,
        print_every_nth_iteration=5,
        use_Sk=False,
        debug=False,
        matrix_debug=False,
    ):
        self.input_dim = input_dim
        # since each surrogate model f_k represents probability p_k, the bounds are 0 to 1
        self.bounds = [(0, 1)] * self.input_dim
        self.min = [(0.) * self.input_dim]
        self.fmin = 0
        self.sd = standard_deviation

        self.with_qutip = with_qutip
        self.use_Sk = use_Sk
        self.print_every_nth_iteration = print_every_nth_iteration
        self.debug = debug
        self.matrix_debug = matrix_debug

    # def Fidelity(self, input_qubits, with_qutip, use_Sk, debug, iteration_counter, angles):
    def Fidelity(self, input_qubits, iteration_counter, angles):
        ''' A single surrogate model for the fidelity '''
        # run the GHZ circuit with the input angles to get the 8x8 pure state density matrix, rho
        if self.with_qutip:
            print('USING QUTIP CIRCUIT')
        rho = density_matrix_from_unitary(
            GHZ_circuit_unitary(
                angles,
                matrix_debug=self.matrix_debug,
            ),
            input_qubits,
            with_qutip=self.with_qutip,
        )
        if self.debug:
            print('Density matrix rho =\n', around(rho, 2))
        Fidelity_func = Fidelity_Sk if self.use_Sk else Fidelity_Pk
        fval = Fidelity_func(rho)

        iteration_count = next(iteration_counter)
        # only print the first and then every % nth iteration
        if iteration_count == 1 or not iteration_count % self.print_every_nth_iteration:
            print(
                '\n' + '-' * 50 + f' iteration #{iteration_count}'
                f'\nGHZ Fidelity: angles = {around(angles[0], 3)}'
                f'\nGHZ Fidelity: angles (in degrees) = {[around(angle * (180 / pi), 1) for angle in angles[0]]}'
                f'\nGHZ Fidelity: fval = f(angles) = {around(fval, 3)}'
            )

        return fval

    def Infidelity(self, *args):
        I = 1 - self.Fidelity(*args)
        print('Infidelity is', round(I, 3))
        return I

    def f_k(self, k, input_qubits, iteration_counter, angles):
        ''' Surrogate model '''
        # NOTE: initial angles specified by `X` kwarg of `BayesianOptimization` class instantiated
        # in running of the script

        # run the GHZ circuit with the input angles to get the 8x8 pure state density matrix, rho
        if self.with_qutip:
            print('USING QUTIP CIRCUIT')
        rho = density_matrix_from_unitary(
            GHZ_circuit_unitary(angles),
            input_qubits,
            with_qutip=self.with_qutip,
        )

        print('rho =\n', rho)

        # evaluate the surrogate model
        measurement_operator = S_operator if self.use_Sk else P_operator
        print(f'USING OPERATOR = {measurement_operator}')
        fval = trace(rho @ measurement_operator(k))

        iteration_count = next(iteration_counter)
        # only print the first and then every % nth iteration
        if iteration_count == 1 or not iteration_count % self.print_every_nth_iteration:
            print(
                '\n' + '-' * 50 + f' iteration #{iteration_count}'
                f'\nGHZ f_{k}: angles = {around(angles[0], 3)}'
                f'\nGHZ f_{k}: angles (in degrees) = {[around(angle * (180 / pi), 1) for angle in angles[0]]}'
                f'\nGHZ f_{k}: fval = f(angles) = {around(fval, 3)}'
            )

        # TODO: do I want to add noise?
        return fval


if __name__ == '__main__':
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

    with_qutip = args.with_qutip  # default False

    random_seed = args.random_seed
    if random_seed != 0:
        random.seed(random_seed)

    print_every_nth_iteration = args.print_every_nth_iteration
    use_Sk = args.use_Sk
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
            angles = [list(random.uniform(0, limit[units], num_angles))]
        elif ghz:
            angles = [list(random.uniform(0, limit[units], num_angles))]
    else:
        angles = [float(angle) for angle in angles]
    angles_pi, angles_radians, angles_degrees = convert_angles(angles, units)

    ghz = GhzStateFidelity(
        with_qutip=with_qutip,
        print_every_nth_iteration=args.print_every_nth_iteration,
        use_Sk=use_Sk,
        debug=debug,
        matrix_debug=matrix_debug,
    )

    func_to_optimize = ghz.Fidelity
    func_to_optimize = ghz.Infidelity


    bopt_kwargs = {
        'domain': [{
            'name': 'angles',
            'type': 'continuous',
            'domain': (0, 2*pi),
            'dimensionality': 6,
        }],

        # this is the default but emphasise all 6 input angles are independent
        # constraints imply inequalities relating the different inputs
        # constraints=None,

        # :model_type: type of model to use as surrogate:
        #     - 'GP', standard Gaussian process.
        #     - 'GP_MCMC', Gaussian process with prior in the hyper-parameters.
        #     - 'sparseGP', sparse Gaussian process.
        #     - 'warperdGP', warped Gaussian process.
        #     - 'InputWarpedGP', input warped Gaussian process
        #     - 'RF', random forest (scikit-learn).
        'model_type': 'GP',

        'X': array(angles_radians),
        # TODO: not sure whether to specify a Y too?
        # Y = func_to_optimize(args.input_qubits, count(1), angles),

        # defaults to 5, is used in generating input angles
        # 'initial_design_numdata': 100,

        # defaults to 'random' - can also be 'grid', 'latin' or 'sobol' - see GPyOpt_designs.py
        'initial_design_type': 'random',


        # TODO: figure out if this is the right acquisition
        # :acquisition_type: type of acquisition function to use.
        #     - 'EI', expected improvement.
        #     - 'EI_MCMC', integrated expected improvement (requires GP_MCMC model).
        #     - 'MPI', maximum probability of improvement.
        #     - 'MPI_MCMC', maximum probability of improvement (requires GP_MCMC model).
        #     - 'LCB', GP-Lower confidence bound.
        #     - 'LCB_MCMC', integrated GP-Lower confidence bound (requires GP_MCMC model).
        'acquisiton_type': 'EI',

        # TODO: investigate DIRECT and CMA-ES algorithms
        # :acquisition_optimizer_type: type of acquisition function to use.
        #     - 'lbfgs', L-BFGS.
        #     - 'DIRECT', Dividing Rectangles.
        #     - 'CMA', covariance matrix adaptation.
        'acquisition_optimizer_type': 'lbfgs',

        # TODO: figure out if I want exact_feval (default is False)
        'exact_feval': True,

        # :param evaluator_type: determines the way the objective is evaluated (all methods are equivalent if the batch size is one)
        #     - 'sequential', sequential evaluations.
        #     - 'random', synchronous batch that selects the first element as in a sequential policy and the rest randomly.
        #     - 'local_penalization', batch method proposed in (Gonzalez et al. 2016).
        #     - 'thompson_sampling', batch method using Thompson sampling.
        'evaluator_type': 'sequential',

        'verbosity': True,
        'verbosity_model': True,

        # NOTE: kernel choice by default seems to be GPy.kern.Matern52 as desired
        # this is the R_2(x) kernel of Eq.22 as applied to Eq.19 (of the paper)
    }
    print(f'Bayesian Optimization kwargs:\n{pformat(bopt_kwargs)}')
    if not args.separate:
        # TRY FIDELITY AS A SINGLE SURROGATE MODEL
        bopt = GPyOpt.methods.BayesianOptimization(
            f=partial(
                func_to_optimize,
                args.input_qubits,
                count(1),
            ),
            maximize=func_to_optimize == ghz.Fidelity,
            **bopt_kwargs,
        )
        bopt.run_optimization(
            max_iter=args.max_iter,
            verbosity=True,
        )
        print('\n\n' + '=' * 100)
        param_names = bopt.model.get_model_parameters_names()
        params = bopt.model.get_model_parameters()
        optimal_angles_degrees = [round(angle * (180 / pi), 1) for angle in bopt.x_opt]
        print(
            f'Optimal angles (in degrees) = {around(optimal_angles_degrees, 3)}'
            f'\nOptimal angles (in radians) = {around(bopt.x_opt, 2)}'
            f'\nOptimal {func_to_optimize.__name__}(angles): {around(bopt.fx_opt, 3)}'
            # this just shows final 6 angles?
            f'\nFinal self.X =\n{around(bopt.X, 2)}'
            f'\nmodel params names = {param_names}, len = {len(param_names)}'
            f'\nmodel params = {around(params)}, len = {len(params[0])}'
        )

        bopt._print_convergence()
        bopt.plot_convergence(filename=f'/qcbo/Infidelity_output_convergence.png')
        bopt.save_report(report_file=f'/qcbo/Infidelity_output_report.txt')

        if debug:
            print('dir bopt', pformat(dir(bopt)))
            print(
                f'what are these?'
                f'\n{bopt.acquisition_type=}'
                f'\n{bopt.acquisition=}'
                f'\n{bopt.constraints=}'
                f'\n{bopt.context=}'
                f'\n{bopt.cost=}'
                f'\n{bopt.cum_time=}'
                f'\n{bopt.evaluator_type=}'
                f'\n{bopt.evaluator=}'
                f'\n{bopt.domain=}'
                f'\n{bopt.model_type=}'
                f'\n{bopt.model=}'
                f'\n{bopt.objective_name=}'
                f'\n{bopt.objective=}'
                f'\n{bopt.space=}'
                f'\n{bopt.suggest_next_locations=}'
                f'\n{bopt.suggested_sample=}'
            )

            param_names = bopt.model.get_model_parameters_names()
            params = bopt.model.get_model_parameters()
            print('dir bopt.model', pformat(dir(bopt.model)))
            print(
                f'what are these model attrs?'
                f'\n{bopt.model.kernel=}'
                f'\n{bopt.model.model=}'
                f'\n{bopt.model.ARD=}'
                f'\n{bopt.model.kernel=}'
                f'\nf min = {bopt.model.get_fmin()}'
                f'\nmodel params names = {param_names}, len = {len(param_names)}'
                f'\nmodel params = {params}, len = {len(params[0])}'
            )
