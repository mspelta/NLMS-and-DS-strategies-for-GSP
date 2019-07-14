% ======================================================================= %
% == COPPE/UFRJ - Programa de Engenharia Eletrica (PEE) ================= %
% == Script: script01_Comparison ======================================== %
% == Responsible: Marcelo Jorge Mendes Spelta - Date: 2019/06/20 ======== %
% == E-mail: marcelo.spelta@smt.ufrj.br ================================= %
% ======================================================================= %

clc;clear;

% ======================================================================= %
% -- Loading original graph signal, bandlimited graph signal, U_f and  -- %
% -- D_s. Graph signal is obtained from the NOAA dataset. --------------- %
load('General_Temperature_Data')
% ======================================================================= %

%%  ==================================================================== %%
% ======================================================================= %
% == I - Comparison among the LMS, RLS and NLMS Algorithm =============== %

% ----------------------------------------------------------------------- %
% -- Initializing mean measurements variables --------------------------- %
mean_error_mat_comp = [];
mean_MSD_mat_comp = [];
mean_elapsedTime_mat_comp = [];

% ----------------------------------------------------------------------- %
% -- Basic Simulation Settings ------------------------------------------ %
numberSamples = 95;  
ensemble = 250;
numbCheckedUpdates = 10; 

% alg_param_vec = [ mu_L beta_R mu_N / mu_L beta_R mu_N / mu_L beta_R mu_N ]
alg_param_vec = [ 1.5 0.5 0.5 ]; 
% alg_selection_vec = [ LMS RLS NLMS / LMS RLS NLMS / LMS RLS NLMS ]
alg_selection_vec = [ 1 2 3 ];  
DS_strategy_vec = ones(1,3); % Not important since we use kappa=0
kappa_vec = zeros(1,3);      % kappa=0 so we don't use the DS strategies
% ----------------------------------------------------------------------- %
% -- Generating diagonal variance vectors. That will be used to form C_w % 
rng(2);
variance_vector_matrix = 0.01*ones(205,3); 
% ----------------------------------------------------------------------- %
% ======================================================================= %

%%  ==================================================================== %%
% ======================================================================= %
% == Running Simulation and Obtaining Results =========================== %

for int_counter = 1:3  
    int_counter
    % ------------------------------------------------------------------- %
    % ------------------------------------------------------------------- %
    [ mean_error_mat_comp(:,int_counter), mean_MSD_mat_comp(:,int_counter), ...
        mean_x_vec_matrix(:,:,int_counter), ~, mean_elapsedTime_mat_comp(:,int_counter) ] = ...
      graph_signal_reconstruction( completeDatasetMatrix(:,1)*ones(1,95) , ... 
            D_s, U_f, diag(variance_vector_matrix(:,int_counter)), kappa_vec(int_counter), ...
            alg_param_vec(int_counter), alg_selection_vec(int_counter), DS_strategy_vec(int_counter), ...
            numberSamples, numbCheckedUpdates, ensemble ) ;
    % ------------------------------------------------------------------- %
    % ------------------------------------------------------------------- %
end

% ======================================================================= %

%%  ==================================================================== %%
% ======================================================================= %
% == Saving important data for plotting results in another script ======= %

save('./Simu_Results/01-simu-StaticComparison_LMS-RLS-NLMS','mean_error_mat_comp','mean_MSD_mat_comp','mean_x_vec_matrix',...
    'mean_elapsedTime_mat_comp','alg_param_vec','alg_selection_vec','ensemble');

% ======================================================================= %
% ======================================================================= %
%%  ==================================================================== %%
% ======================================================================= %
% == Running Simulation and Obtaining Results =========================== %

for int_counter = 1:3 
    int_counter
    % ------------------------------------------------------------------- %
    % ------------------------------------------------------------------- %
    [ mean_error_mat_comp(:,int_counter), mean_MSD_mat_comp(:,int_counter), ...
        mean_x_vec_matrix(:,:,int_counter), ~, mean_elapsedTime_mat_comp(:,int_counter) ] = ...
      graph_signal_reconstruction( completeDatasetMatrix , ... 
            D_s, U_f, diag(variance_vector_matrix(:,int_counter)), kappa_vec(int_counter), ...
            alg_param_vec(int_counter), alg_selection_vec(int_counter), DS_strategy_vec(int_counter), ...
            numberSamples, numbCheckedUpdates, ensemble ) ;
    % ------------------------------------------------------------------- %
    % ------------------------------------------------------------------- %
end

% ======================================================================= %

%%  ==================================================================== %%
% ======================================================================= %
% == Saving important data for plotting results in another script ======= %

save('./Simu_Results/01-simu-DynamicComparison_LMS-RLS-NLMS','mean_error_mat_comp','mean_MSD_mat_comp','mean_x_vec_matrix',...
    'mean_elapsedTime_mat_comp','alg_param_vec','alg_selection_vec','ensemble');

% == END OF SCRIPT ====================================================== %
% ======================================================================= %