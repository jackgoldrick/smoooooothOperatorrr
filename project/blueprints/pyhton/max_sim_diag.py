from p_power import *
import numpy as np
import torch as tc

def max_sim_diag(
    diags_A, matrix_U, p, step=0.1, s_max_diag=50, 
    s_max_power=100, err_a=1e-6, err_scale=100,
):
    # dim2 is number of blocks. each block is a dim1 x dim1 matrix
    dim2, dim1 = diags_A.shape
    diags_A = diags_A.unsqueeze(0).unsqueeze(-1).to(device)
    matrix_U = matrix_U.to(device)
    matrix_U_inv = tc.linalg.inv(matrix_U)
    # during iterations, dims are (s_max_diag, dim2, dim1, dim1)
    step = step * tc.ones(s_max_diag, 1, 1, 1).to(device)
    v_max = tc.zeros(s_max_diag, dim1 * dim2, 1).to(device)
    old_guess = 0
    q = 1 / (1 - 1/p)
    
    stack_A = matrix_U @ (diags_A * matrix_U_inv)
    stack_A = stack_A.flatten(1, 2)
    norm_stack_A, _ = p_power(stack_A, p, err_a=err_a, s_max=s_max_power)
    
    diags_B = tc.randn(s_max_diag, dim2, dim1, 1, dtype=tc.cfloat).to(device)
    w = tc.randn(s_max_diag, 1, dim1, 1, dtype=tc.cfloat).to(device)
    
    converged = tc.zeros(s_max_diag, 1, 1, 1).to(device).bool()
    norm_stack_B = tc.zeros(s_max_diag, 1, 1, 1).to(device)
    v_max = tc.zeros(s_max_diag, dim1 * dim2, 1).to(device)
     
    while True:
        stack_B = matrix_U @ (diags_B * matrix_U_inv)
        stack_B = stack_B.transpose(1, 2).flatten(2, 3)
        
        flat_converged = converged.flatten()
        norm_stack_B[~flat_converged,:,0,0], v_max_tmp = p_power(
            stack_B[~flat_converged,:,:], p,
            v_init=v_max[~flat_converged,:,:],
            err_a=err_a, s_max=s_max_power
        )
        if v_max_tmp is not None:
            v_max[~flat_converged,:,:] = v_max_tmp
        
        diags_B /= norm_stack_B
        diags_sum = tc.sum(diags_A * diags_B, dim=1, keepdim=True)
        v = dual(
            (((w.mH @ matrix_U) * diags_sum.mT) @ matrix_U_inv).mH, 
            q, dim=2
        )
        # import pdb; pdb.set_trace()
        w = dual( matrix_U @ (diags_sum * (matrix_U_inv @ v)), p, dim=2)
        wH_U = w.mH @ matrix_U
        Uinv_v = matrix_U_inv @ v
        
        guess = wH_U @ (diags_sum * Uinv_v)
        guess = tc.abs(guess)
        
        max_guess = tc.max(guess)
        
        if max_guess > norm_stack_A + err_scale * err_a:
            print("Matrix Nrom inequality violated")
            break
        
        converged |= (abs(guess - norm_stack_A) <= err_a) | (guess - old_guess <= err_a)
        
        if tc.all(converged):
            break
        
        gradient = wH_U.mT * diags_A * Uinv_v
        
        step[converged] = 0
        diags_B += step * gradient
        
        old_guess = guess
            
    if tc.min(norm_stack_A - max_guess) < -err_scale * err_a:
        print("Matrix norm inequality violated")
        
    return max_guess, norm_stack_A
            