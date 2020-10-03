# Quantum Control with Bayesian Optimization
#### qcbo - said like 'placebo' but 'qucebo'?

Extending work done in https://arxiv.org/pdf/1909.01229.pdf
___
### Setup
1. Install `docker`: https://docs.docker.com/get-docker/
2. Clone this repo
3. Run `make qcbo-dev`

You should now be in a Docker container with `qutip` and `GPyOpt` code setup for development.
(This means you can edit the forked submodules of both of those repos in your local git clone and your scripts will use those code changes :D)

### Run
Run the script `python GHZ_state/GHZ_bayesian_optimization.py` to see the Bayesian optimization process at work in calculating the GHZ Fidelity outlined in Eq.27 of the paper linked above.
