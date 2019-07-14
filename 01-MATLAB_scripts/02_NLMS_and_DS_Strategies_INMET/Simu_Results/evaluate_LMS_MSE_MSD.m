function [MSE, MSD, MSD_Lor] = evaluate_LMS_MSE_MSD(mu, D, U_f, C_w)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    
    A = U_f'*D*U_f;
    B = U_f'*D*C_w*D*U_f;
    
    S_n = zeros(size(A));
    
    for counter = 1:10
        C = mu*(B + A*S_n*A);
        S_n = sylvester(A,A,C);
    end
    
    for counter = 1:299
        if (D(counter,counter) == 1)
            m = U_f(counter,:)';
            est_var_vec(counter,1) = C_w(counter,counter) +  m' * S_n * m;
        else
            est_var_vec(counter,1) = 0;
        end
    end

    %figure;stem(est_var_vec)

    MSD_Lor = (mu/2)*trace( inv(A)*B );
    MSD = trace(S_n);
    MSE = sum(est_var_vec);

    % --------------------------------- 

%     inv_mat = inv(U_f'*D*U_f);
% 
%     for counter = 1:size(C_w,1)
%         if (D(counter,counter) == 1)
%             m = U_f(counter,:)';
%             est_var_vec(counter,1) = C_w(counter,counter) + mu / ( 2 - mu) * m' * inv_mat * (U_f'*D*C_w*D*U_f) * inv_mat*  m;
%         else
%             est_var_vec(counter,1) = 0;
%         end
%     end
% 
%     MSD = mu/(2-mu)* trace( inv_mat * (U_f'*D*C_w*D*U_f) * inv_mat );
%     MSE = sum(est_var_vec);

end

