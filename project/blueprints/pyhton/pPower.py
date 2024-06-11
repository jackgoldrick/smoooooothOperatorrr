import torch as tc
import numpy as np
device = tc.device('cuda' if tc.cuda.is_available() else 'cpu')
tc.set_grad_enabled(False)

def dual(y, p, dim):
    '''
    Returns the Holder p-dual of y along dimension dim
    '''
    if p == 1:
        return tc.sgn(y)
    elif np.isinf(p) or np.isnan(p):
        y_dual = tc.zeros(y.shape).to(device)
        _, ind = y.abs().max(dim=dim, keepdim=True)
        range_shape = [1 for _ in range(len(y.shape))]
        range_shape[dim] = -1
        mask = tc.arange(y.shape[dim]).reshape(range_shape) == ind
        y_dual[mask] = tc.sgn(y[mask])
        return y_dual
    else:
        q = 1 / (1 - 1/p)
        y_dual = (tc.abs(y) ** (p-1) ) * tc.sgn(y)
        return y_dual / tc.linalg.vector_norm(y, ord=p, dim=dim, keepdim=True)
          

def pPower():
    pass

