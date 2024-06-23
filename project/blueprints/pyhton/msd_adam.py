from p_power import *
import numpy as np
import torch as tc

def max_sim_diag(
    diags_A, matrix_U, p, step_size=0.1, s_max_diag=50, 
    s_max_power=100, err_a=1e-7, err_scale=100, lyfe_cycle=0, decay=[0.999, 0.9], max_breaders=4
):
    assert 0 < max_breaders < s_max_diag, "max_breaders must bee less than _s_max_diag"
    violations = 0
    generation = 1
    og_gen = True
    max_guess = tc.Tensor(lyfe_cycle+1, 1).to(device)
    # dim2 is number of blocks. each block is a dim1 x dim1 matrix
    dim2, dim1 = diags_A.shape
    diags_A = diags_A.unsqueeze(0).unsqueeze(-1).to(device)
    matrix_U = matrix_U.to(device)
    matrix_U_inv = tc.linalg.inv(matrix_U)
    # during iterations, dims are (s_max_diag, dim2, dim1, dim1)
    step = step_size * tc.ones(s_max_diag, 1, 1, 1).to(device)
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
        new_birth = False
        step_adj = step * (np.sqrt(1 - decay[0]**t) / (1 - decay[1]**t))
        
    
        stack_B = matrix_U @ (diags_B * matrix_U_inv)
        stack_B = stack_B.transpose(1, 2).flatten(2, 3)
        
        # import pdb; pdb.set_trace()
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
            
        # else:
        #     stacks = (stack_A @ stack_B)
        #     v = dual(w @ stacks, q, dim=2)
        #     w = dual(stacks @ v, p, dim=2)
        #     guess = w @ stacks @ v
            
        #     del stacks  # Delete the tensor reference

        #     # Clear cache if the tensor was on GPU
        #     if tc.cuda.is_available():
        #         tc.cuda.empty_cache()
            
        #     Uinv_v = matrix_U_inv @ v
        
        guess = tc.abs(guess)
        
        max_guess[generation-1] = tc.max(guess)
        
        
        
            # direction *= -1
            
        
        converged |= (abs(guess - norm_stack_A) <= err_a) | ((tc.abs(guess - old_guess)) <= err_a * err_scale)
        
        if tc.all(converged):
            if (max_guess[generation-1, 0] > norm_stack_A + err_scale * err_a):
                violations += 1
            
            if generation > lyfe_cycle:
                
                break
            res = tc.topk(guess[:,0,0,0], max_breaders)
            locations = res.indices
            
            parents = diags_B[locations, :, :, :]

            # import pdb; pdb.set_trace()
            # mitosis = (parents * parents.mT)
            # diags_B = (mitosis[:,:,:,1] * mitosis[:,:,:,0]) ^ 1/p 
            # 
            # print(diags_B.shape)
            diags_B = (parents * parents.transpose(1,2))
            
            # import pdb; pdb.set_trace()
            generation += 1
            t = 1
            # converged = (0 * converged).bool() 
            converged = tc.zeros(max_breaders, 1, 1, 1, dtype=tc.cfloat).to(device).bool()
            # import pdb; pdb.set_trace()
            old_guess = (res.values).unsqueeze(-1).unsqueeze(-1).unsqueeze(-1)  
            
            if og_gen:
                del v_max, norm_stack_B, w, step,
                A = A[max_breaders-1, :, :, :]
                F = F[max_breaders-1, :, :, :]
                
                v_max = tc.randn(max_breaders, dim1 * dim2, 1, dtype=tc.cfloat).to(device)
                norm_stack_B = tc.ones(max_breaders, 1, 1, 1).to(device)
                w = tc.randn(max_breaders, 1, dim1, 1, dtype=tc.cfloat).to(device)
                step = step_size * tc.ones(max_breaders, 1, 1, 1).to(device)
                og_gen = False
            new_birth = True
            
            

        if (not new_birth):
            gradient = wH_U.mT * diags_A * Uinv_v
            step[converged] = 0
            
            
            
            
            A = decay[0] * A + (1 - decay[0]) * ((gradient ** 2))
            
            F = decay[1] * F + (1 - decay[1]) * (gradient)
            
            
            diags_B += (((step_adj) / (tc.sqrt(A) + err_a  / err_scale) ) * F)
            
            t += 1
            old_guess = guess
            
    final_max_guess = tc.max(max_guess)
    if tc.min(norm_stack_A - final_max_guess) < -err_scale * err_a:
        
        print("Value returned will violate the Matrix norm inequality")
        print(f'Number of violations for p={p:.2f}: {violations}')
        
    return final_max_guess, norm_stack_A
            