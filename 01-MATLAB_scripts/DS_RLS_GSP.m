% ======================================================================= %
% -- Responsible: Marcelo J. M. Spelta - Date: 2018/02/22
% -- Algorithm 1 Function
% -- Based on the LMS algorithm for Graph Signals presented in the paper
% -- called "Adaptive LMS (Least Mean Squares) Estimation of Graph Signals"
% ======================================================================= %

% ======================================================================= %

function [ new_x_hat, PSI, update_flag, elapsedTime ] = ...
    DS_RLS_GSP( x_w, x_hat, U_f, D, M_prime, B_r, beta, previous_PSI, var_vec, DS_strategy, kappa )
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
      PSI = beta*previous_PSI + M_prime;
      %new_x_hat = x_hat + ( U_f * ( inv(PSI) * ( B_r * error_vec ) ) );  % %new_x_hat = x_hat + U_f * ( PSI \ ( U_f' * D * ( C_w \ error_vec) ) ) ;
      new_x_hat = x_hat + U_f * ( PSI \ ( B_r * error_vec) ) ;
      elapsedTime = toc;
  else
      PSI = previous_PSI;
      new_x_hat = x_hat;
      elapsedTime = 0;
  end
  % ===================================================================== %
    
end

% ======================================================================= %

