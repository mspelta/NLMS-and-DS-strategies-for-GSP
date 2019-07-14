function [MSE, MSD] = evaluate_NLMS_MSE_MSD(mu, D, U_f, C_w, x_diff)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    inv_mat = inv(U_f'*D*U_f);

    for counter = 1:size(C_w,1)
        if (D(counter,counter) == 1)
            m = U_f(counter,:)';
            est_var_vec(counter,1) = C_w(counter,counter) + mu / ( 2 - mu) * m' * inv_mat * (U_f'*D*C_w*D*U_f) * inv_mat*  m;
        else
            est_var_vec(counter,1) = 0;
        end
    end

    %MSD = mu/(2-mu)* trace( inv_mat * (U_f'*D*C_w*D*U_f) * inv_mat );
    P = D*U_f*inv(U_f'*D*U_f);
    MSD = x_diff'*P*P'*x_diff + x_diff'*x_diff - 2*x_diff'*U_f*inv(U_f'*D*U_f)*U_f'*D*x_diff + (mu/(2-mu))*trace(P'*C_w*P);
    %MSE = sum(est_var_vec);
    MSE = x_diff'*D*x_diff + trace(D*C_w) - 2*x_diff'*D*U_f*inv(U_f'*D*U_f)*U_f'*D*x_diff + ...
        trace(U_f'*D*(x_diff*x_diff' + mu*C_w/(2-mu))*D*U_f*inv(U_f'*D*U_f));
    
end

