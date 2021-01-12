#!/usr/bin python3
'''
Different types of initial_design
https://nbviewer.ipython.org/github/SheffieldML/GPyOpt/blob/master/manual/GPyOpt_initial_design.ipynb
'''

import numpy as np
import GPyOpt
import GPy
from GPyOpt.experiment_design import initial_design
import matplotlib.pyplot as plt


func  = GPyOpt.objective_examples.experimentsNd.alpine1(input_dim=2)

mixed_domain =[{'name': 'var1_2', 'type': 'continuous', 'domain': (-10,10),'dimensionality': 1},
            {'name': 'var5', 'type': 'continuous', 'domain': (-1,5)}]

space = GPyOpt.Design_space(mixed_domain)
data_init = 500

if __name__ == '__main__':

    for design_name in ('random', 'grid', 'sobol', 'latin'):
        X = initial_design(design_name, space, data_init)
        plt.plot(X[:,0],X[:,1],'b.')
        plt.savefig(f'/qcbo/design_{design_name}.png')
        plt.close()
