function [ update_flag ] = update_checking_test( error_vec, S, var_vec, DS_strategy, kappa )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    update_flag = 0;  % initializing update_flag variable -> No update
    norm_error_vec = zeros(size(error_vec));
    
    % =================================================================== %
    % == UPDATE CHECKING PART =========================================== %
    if(DS_strategy == 1) % Component-wise error inspection
        % --------------------------------------------------
        ref_gamma_vec = kappa * sqrt(var_vec);
        % --------------------------------------------------
        for i = 1:size(error_vec)
            if( abs( error_vec(i) ) > ref_gamma_vec(i) )
                update_flag = 1;
            end
        end
        % --------------------------------------------------
    else                 % Norm l-2 error inspection
        % --------------------------------------------------
        for i = 1:size(error_vec)
            if( var_vec(i) ~= 0 )
                norm_error_vec(i,1) = error_vec(i) / (sqrt(var_vec(i)));
            else
                norm_error_vec(i,1) = 0;
            end
        end
        % --------------------------------------------------
        if( norm_error_vec'*norm_error_vec > kappa*S )   % 0.5 * kappa * 210
                update_flag = 1;
        end
        % --------------------------------------------------
    end

end

