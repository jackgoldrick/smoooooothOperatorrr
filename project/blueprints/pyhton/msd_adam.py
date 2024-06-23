from p_power import *
import numpy as np
import torch as tc

def max_sim_diag(
    diags_A, matrix_U, p, step=0.1, s_max_diag=50, 
    s_max_power=100, err_a=1e-6, err_scale=100, lyfe_cycle=0, decay=[0.999, 0.9], max_breaders=4
):
    assert 0 < max_breaders < s_max_diag, "max_breaders must bee less than _s_max_diag"
    violations = 0
    generation = 1
    # dim2 is number of blocks. each block is a dim1 x dim1 matrix
    dim2, dim1 = diags_A.shape
    diags_A = diags_A.unsqueeze(0).unsqueeze(-1).to(device)
    matrix_U = matrix_U.to(device)
    matrix_U_inv = tc.linalg.inv(matrix_U)
    # during iterations, dims are (s_max_diag, dim2, dim1, dim1)
    step = step * tc.ones(s_max_diag, 1, 1, 1).to(device)
    v_max = tc.zeros(s_max_diag, dim1 * dim2, 1).to(device)
    old_guess = 0
    q = h_conj(p)
    
    stack_A = matrix_U @ (diags_A * matrix_U_inv)
    stack_A = stack_A.flatten(1, 2)
    norm_stack_A, _ = p_power(stack_A, p, err_a=err_a, s_max=s_max_power)
    # print(norm_stack_A)
    diags_B = tc.randn(s_max_diag, dim2, dim1, 1, dtype=tc.cfloat).to(device)
    w = tc.randn(s_max_diag, 1, dim1, 1, dtype=tc.cfloat).to(device)
    
    converged = tc.zeros(s_max_diag, 1, 1, 1, dtype=tc.cfloat).to(device).bool()
    v_max = tc.randn(s_max_diag, dim1 * dim2, 1, dtype=tc.cfloat).to(device)
    norm_stack_B = tc.ones(s_max_diag, 1, 1, 1).to(device)


    t = 1
    A = tc.zeros(s_max_diag, 1, 1, 1).to(device)
    F = tc.zeros(s_max_diag, 1, 1, 1).to(device)

    while True:
        step_adj = step * (np.sqrt(1 - decay[0]**t) / (1 - decay[1]**t))
        
        if ~lyfe_cycle:
            stack_B = matrix_U @ (diags_B * matrix_U_inv)
            stack_B = stack_B.transpose(1, 2).flatten(2, 3)
       
        
        flat_converged = converged.flatten()
        tc.nn.init.ones_(norm_stack_B)
        norm_stack_B[~flat_converged,0,:,:], v_max_tmp = p_power(
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
        w = dual( matrix_U @ (diags_sum * (matrix_U_inv @ v)), p, dim=2)
        wH_U = w.mH @ matrix_U
        Uinv_v = matrix_U_inv @ v
        
        guess = wH_U @ (diags_sum * Uinv_v)
        guess = tc.abs(guess)
        
        max_guess = tc.max(guess)
        
        
        if max_guess > norm_stack_A + err_scale * err_a:
            violations += 1
            # direction *= -1
            
        
        converged |= (abs(guess - norm_stack_A) <= err_a) | ((tc.abs(guess - old_guess)) <= err_a * err_scale)
        
        if tc.all(converged):
            if generation > lyfe_cycle:
                
                break
            res = tc.topk(guess[:,0,0,0], max_breaders)
            locations = res.indices
            
            parents = stack_B[locations, :, :]

            # import pdb; pdb.set_trace()
            children = (parents * parents.transpose(0,2))
            
            # import pdb; pdb.set_trace()
            generation += 1
            stack_B = children
            # converged = (0 * converged).bool() 
            converged = tc.zeros(s_max_diag, 1, 1, 1, dtype=tc.cfloat).to(device).bool()
            

        if ~lyfe_cycle:
            gradient = wH_U.mT * diags_A * Uinv_v
            step[converged] = 0
            
            
            
            
            A = decay[0] * A + (1 - decay[0]) * ((gradient ** 2))
            
            F = decay[1] * F + (1 - decay[1]) * (gradient)
            
            
            diags_B += (((step_adj) / (tc.sqrt(A) + err_a  / err_scale) ) * F)
            
            
        
        
            t += 1
        old_guess = guess 
    if tc.min(norm_stack_A - max_guess) < -err_scale * err_a:
        
        print("Value returned will violate the Matrix norm inequality")
        print(f'Number of violations for p={p:.2f}: {violations}')
        
    return max_guess, norm_stack_A
            