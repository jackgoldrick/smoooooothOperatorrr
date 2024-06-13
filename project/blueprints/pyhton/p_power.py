import torch as tc
import numpy as np
device = tc.device('cuda' if tc.cuda.is_available() else 'cpu')
tc.set_grad_enabled(False)
print(f'Using device: {device}')
def dual(y, p, dim):
    '''
    Returns the Holder p-dual of y along dimension dim
    '''
    if p == 1:
        return y.sgn()
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
        return y_dual / tc.linalg.vector_norm(y_dual, ord=q, dim=dim, keepdim=True)
          

def p_power(matrix, p, err_a=1e-6, s_max=25, v_init=None):
    '''
    Matrix operator p-norm along last two dimensions of matrix.
    '''
    # We want dims to be (batch, s_max, n, m)
    
    if len(matrix.shape) == 2:
        matrix = matrix.unsqueeze(0)  # add batch dim
    matrix = matrix.unsqueeze(1)  # add s_max dim
    dims = list(matrix.shape)
    dims[1] = s_max
    
    if p == 1 or np.isinf(p):
        return tc.linalg.matrix_norm(matrix, ord=p), None
   
    matrix_prime = matrix.conj().transpose(2, 3)
    
    q = 1 / (1 - 1/p)
    
    x_dims = (dims[0], dims[1], dims[3], 1)
    x = tc.complex(
        *[tc.normal(mean=tc.zeros(x_dims), std=1) for _ in range(2)]
    ).to(device)
    
    if v_init is not None:
        x = tc.concat([x, v_init], dim=1)
        
    x_norm = tc.linalg.vector_norm(x, ord=p, dim=2, keepdims=True)
    x /= x_norm
    old_guess = 0
    
    while(True):
        y = matrix @ x
        y_dual = dual(y, p, dim=2)
        # import pdb; pdb.set_trace()
        guess, ind = tc.max(
            tc.linalg.vector_norm(y, ord=p, dim=2, keepdim=True),
            dim=1,
            keepdim=True,
        )
        # ind = ind.flatten()
        # guess = guess.flatten()
        # print(max(abs(guess - old_guess)))
        if tc.max(abs(guess - old_guess)) < err_a:
            ran = tc.arange(s_max)
            v_max = x[:,ind.flatten() == ran, :, :]
            break
        
        z = matrix_prime @ y_dual
        x = dual(z, q, dim=2)
        old_guess = guess
        
        if guess.isnan().any():
            # import pdb; pdb.set_trace()
            print('Floating point error!')
            return 0, None
        
    return guess, v_max
    
    
    


