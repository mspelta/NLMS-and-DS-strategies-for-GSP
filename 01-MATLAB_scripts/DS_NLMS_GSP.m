% ======================================================================= %
% -- Responsible: Marcelo J. M. Spelta - Date: 2018/02/22
% -- Algorithm 2 Function
% -- Based on the algorithm obtained by minimizing the norm ||w(k+1) - 
% -- w(k)||_2^2 considering the constraints y(k) - DBw(k+1) = 0 and
% -- Bw(k+1) = w(k+1)
% ======================================================================= %

% ======================================================================= %
function [ new_x_hat, update_flag, elapsedTime ] = ...
    DS_NLMS_GSP( x_w, x_hat, D, mu_Bn_matrix, var_vec, DS_strategy, kappa )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    error_vec = D*(x_w - x_hat);
    
    % =================================================================== %
    % == UPDATE CHECKING PART =========================================== %
    update_flag = update_checking_test( error_vec, sum(diag(D)), var_vec, DS_strategy, kappa );
    
    % ===================================================================== %
    % == TRADITIONAL NLMS ALGORITHM ======================================= %
    if ( update_flag == 1 )  
        tic;
        new_x_hat = x_hat + mu_Bn_matrix * error_vec ;
        elapsedTime = toc;
    else
        new_x_hat = x_hat;
        elapsedTime = 0;
    end
    % ===================================================================== %
end
% ======================================================================= %