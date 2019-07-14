function [MSE, MSD] = evaluate_RLS_MSE_MSD(beta, D, U_f, C_w)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    
    inv_mat = inv(U_f'*D*inv(C_w)*D*U_f);

    S_r = inv_mat * ( U_f'*D*inv(C_w)*D*C_w*D*inv(C_w)*U_f ) * inv_mat;

    for counter = 1:size(C_w,1)
        if (D(counter,counter) == 1)
            m = U_f(counter,:)';
            est_var_vec(counter,1) = C_w(counter,counter) + (1-beta)/(1+beta) * m' * S_r * m;
        else
            est_var_vec(counter,1) = 0;
        end
    end
               
    
    MSD = (1-beta)/(1+beta) * trace( S_r );
    MSE = sum(est_var_vec);


end

