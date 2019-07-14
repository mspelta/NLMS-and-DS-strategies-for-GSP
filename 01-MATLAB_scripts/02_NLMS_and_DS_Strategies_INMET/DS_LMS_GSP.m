% ======================================================================= %
% -- Responsible: Marcelo J. M. Spelta - Date: 2018/02/22
% -- Algorithm 1 Function
% -- Based on the LMS algorithm for Graph Signals presented in the paper
% -- called "Adaptive LMS (Least Mean Squares) Estimation of Graph Signals"
% ======================================================================= %

% ======================================================================= %

function [ new_x_hat, update_flag, elapsedTime ] = DS_LMS_GSP( x_w, x_hat, D, mu_Bl_matrix, var_vec, DS_strategy, kappa )
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
        new_x_hat = x_hat + mu_Bl_matrix * error_vec;
        elapsedTime = toc;
    else
        new_x_hat = x_hat;
        elapsedTime = 0;
    end
  % ===================================================================== % 
end

% ======================================================================= %

