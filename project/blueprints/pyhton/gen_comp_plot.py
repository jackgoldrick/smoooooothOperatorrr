%matplotlib widget
import p_power as pp # Before you ask, yes, I know it's a bad name. I'm working on it.... OMG THAT WAS COPILOT'S SUGGESTION! I LOVE COPILOT!
import matplotlib.pyplot as plt
import numpy as np
import torch as tc
import time
from typing import TypeAlias
from scipy.linalg import hadamard
import csv

device = tc.device('cuda' if tc.cuda.is_available() else 'cpu')
tc.set_grad_enabled(False)
print(f'Using device: {device}')


def plot(method='p', kind=None, file_path=None, dims=None, dp=.1, p_max=25, s_max= 25, pq=False):

    if file_path == None:

        if type(kind) == str:
            
            if dims is not None:
                if len(dims) == 1:
                    row_dim, col_dim = dims 
                elif len(dims) == 2:
                    row_dim, col_dim = dims[0], dims[1]
                else:
                    raise ValueError('dims must be a list with one or two elements')
                    
            match type:
            
                case 'r':
                    matrix = tc.randn(row_dim, col_dim, dtype=tc.cfloat).to(device) 

                case 'c':
                    matrix = tc.complex(tc.randn(row_dim, col_dim, dtype=tc.cfloat),  tc.randn(row_dim, col_dim, dtype=tc.cfloat)).to(device)
                    
                case 'h': 
                    matrix = tc.tensor(hadamard(col_dim), dtype= tc.cfloat)    

                case _:
                    raise ValueError('Random Contrustion Type not availabe. If you are looking to add a random contrstruction please contact Jack Goldrick (@jackgoldrick) with request')
                
                
        elif type(kind) == list:
            matrix = tc.tensor(kind, dtype=tc.float)
                
    else:
        matrix = csv


