# Quantum Control with Bayesian Optimization
#### qcbo - said like 'placebo' but 'qucebo'?

Extending work done in https://arxiv.org/pdf/1909.01229.pdf
___
### Setup
1. Install `docker`: https://docs.docker.com/get-docker/
2. Clone this repo
3. Run `make qcbo-dev`

You should now be in a Docker container with `qutip` and `GPyOpt` code setup for development.
(This means you can edit the forked submodules of both of those repos in your local git clone and your scripts will use those code changes 😄).

`qiskit` is also installed, although their is no editable code mounted here. Judging from their
[installation instructions](https://qiskit.org/documentation/contributing_to_qiskit.html#install-install-from-source-label), they're in dire need of a Dockerfile...

### Scripts
#### calculate_fidelity.py
e.g. 1)
`python circuit_calculations/calculate_fidelity.py --angles 0.5 0 0 0 0 0 --with_qiskit --ghz`
This will calculate the fidelity of the state produced by:
```
             ┌─────────┐          ┌───────┐
        q_0: ┤ RX(π/2) ├──■────■──┤ RX(0) ├
             └┬───────┬┘┌─┴─┐  │  ├───────┤
        q_1: ─┤ RX(0) ├─┤ X ├──┼──┤ RX(0) ├
              ├───────┤ └───┘┌─┴─┐├───────┤
        q_2: ─┤ RX(0) ├──────┤ X ├┤ RY(0) ├
              └───────┘      └───┘└───────┘
```

e.g. 2)
`python circuit_calculations/calculate_fidelity.py --angles 0 60 90 120 180 240 270 300 --units degrees --permute --bell`
This will iterate through all 8^4 = 4096 permutations of the 4 angles provided for a quantum circuit comprised of rotation and control gates, returning information on the output state's fidelity to the Bell state.

e.g. a perfect fidelity state is given by
```
        ┌─────────┐      ┌─────────┐
   q_0: ┤ RX(π/2) ├───■──┤ RX(π/2) ├
        ├─────────┴┐┌─┴─┐├─────────┤
   q_1: ┤ RX(3π/2) ├┤ X ├┤ RY(π/2) ├
        └──────────┘└───┘└─────────┘
```

Note: calculations are done by default directly through numpy arrays (a lot quicker) unless otherwise specified as `--with_qutip` (a lot slower) or `--with_qiskit` (in between). Future work will hopefully add a `--with_cirq` option, and maybe even GPU accelerated array calculations (if I can figure out PyOpenCL!).

#### bayesian_optimize.py

Run the script `python circuit_calculations/bayesian_optimize.py` to see the Bayesian optimization process at work in calculating the GHZ Fidelity outlined in Eq.27 of the paper linked above.
