import ipympl
import p_power as pp # Before you ask, yes, I know it's a bad name. I'm working on it.... OMG THAT WAS COPILOT'S SUGGESTION! I LOVE COPILOT!
import matplotlib.pyplot as plt
import numpy as np
import torch as tc
import time
from typing import TypeAlias
from scipy.linalg import hadamard
import pandas as pd

device = tc.device('cuda' if tc.cuda.is_available() else 'cpu')
tc.set_grad_enabled(False)
print(f'Using device: {device}')


def p_plot(kind=None, file_path=None, dims=None, dp=.1, p_max=25, s_max= 25, pq=False):

    if file_path == None:

        if type(kind) == str:
            
            if dims is not None:
                if len(dims) == 1:
                    row_dim, col_dim = dims[0], dims[0]
                elif len(dims) == 2:
                    row_dim, col_dim = dims[0], dims[1]
                else:
                    raise ValueError('dims must be a list with one or two elements')
                    
            match kind:
            
                case 'r':
                    matrix = tc.randn(row_dim, col_dim, dtype=tc.cfloat).to(device) 

                case 'c':
                    matrix = tc.complex(tc.randn(row_dim, col_dim),  tc.randn(row_dim, col_dim)).to(device)
                    
                case 'h': 
                    matrix = tc.tensor(hadamard(col_dim), dtype= tc.cfloat).to(device)

                case _:
                    raise ValueError('Random Contrustion Type not availabe. If you are looking to add a random contrstruction please contact Jack Goldrick (@jackgoldrick) with request')
                
                
        elif type(kind) == list:
            matrix = tc.tensor(kind, dtype=tc.float)
                
    else:
        matrix = tc.tensor(pd.read_csv(file_path, delimiter=','), dtype=tc.cfloat)
        
    p = np.arange(1, p_max, dp)
    norms = np.zeros(p.shape)
    obj_max =None
    for j in p:
        norms, obj_max = pp.p_power(matrix=matrix, p=j, s_max=s_max, v_init=obj_max)
        
        
        
    f1 = plt.figure()
    plt.xlabel('p')
    plt.ylabel('Norm')
    plt.plot(p, norms, label='Power Method')    
    
        
        

