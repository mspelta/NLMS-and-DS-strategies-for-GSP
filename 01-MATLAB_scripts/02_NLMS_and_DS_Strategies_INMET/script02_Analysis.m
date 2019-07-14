% ======================================================================= %
% == COPPE/UFRJ - Programa de Engenharia Eletrica (PEE) ================= %
% == Script: script02_Analysis ========================================== %
% == Responsible: Marcelo Jorge Mendes Spelta - Date: 2019/03/29 ======== %
% == E-mail: marcelo.spelta@smt.ufrj.br ================================= %
% ======================================================================= %

clc;clear;

% ======================================================================= %
% -- Loading original graph signal, bandlimited graph signal, U_f and  -- %
% -- D_s. Graph signal is obtained from the INMET dataset. -------------- %
load('General_Bandlimited_GS_Data')
% ======================================================================= %

%% ===================================================================== %%
% == Initializing general variables used in all simulations ============= %

ensemble = 5; %10; %1000; %1000; 
numbCheckedUpdates = 1000;
graph_signal_type = 2; %1;

DS_strategy_vec = ones(1,8);    % Not important since we use kappa=0
kappa_vec = zeros(1,8);         % kappa=0 so we don't use the DS strategies

% ----------------------------------------------------------------------- %
% -- Generating diagonal variance vectors. That will be used to form C_w % 
variance_vector_matrix = [ 0.01*ones(299,4) ]; 
rng(2);
variance_vector = 0.005 + 0.010*rand(299,1);
variance_vector_matrix = [ variance_vector_matrix variance_vector*ones(1,4) ];
% ----------------------------------------------------------------------- %
% ======================================================================= %

%%  ==================================================================== %%
% ======================================================================= %
% == I - NLMS Algorithm (2019-01-06) ==================================== %

% ----------------------------------------------------------------------- %
% -- Initializing mean measurements variables --------------------------- %
mean_error_mat_LMS_comp = [];
mean_MSD_mat_LMS_comp = [];
mean_orig_MSD_mat_LMS_comp = [];
update_counter_vec_LMS_comp = [];

% ----------------------------------------------------------------------- %
% -- Basic Simulation Settings ------------------------------------------ %
numberSamples = 3000; 

mu_vec = [ 0.05 0.10 0.25 0.50 0.05 0.10 0.25 0.50 ]; % mu_N
alg_selected = 3;                                     % NLMS

% ======================================================================= %

%%  ==================================================================== %%
% ======================================================================= %
% == Running Simulation and Obtaining Results (I) ======================= %

for int_counter = 1:length(mu_vec) 
    int_counter

    [ mean_error_mat_LMS_comp(:,int_counter), mean_MSD_mat_LMS_comp(:,int_counter), ...
      mean_orig_MSD_mat_LMS_comp(:,int_counter), ~, update_counter_vec_LMS_comp(int_counter), ~] = ...
      graph_signal_reconstruction( graph_signal , bandlimited_graph , [ 1 1 ], ...
            D_s, U_f, diag(variance_vector_matrix(:,int_counter)), kappa_vec(int_counter), ...
            mu_vec(int_counter), alg_selected, DS_strategy_vec(int_counter), numberSamples, numbCheckedUpdates, ensemble, graph_signal_type) ;
end

%%  ==================================================================== %%
% ======================================================================= %
% == Saving important data (I) for plotting results in another script === %
save('./Simu_Results/02-simu-DS-NLMS','mean_MSD_mat_LMS_comp',...
    'mean_error_mat_LMS_comp','mu_vec','variance_vector_matrix','ensemble');


%%  ==================================================================== %%
% ======================================================================= %
% == II - LMS Algorithm (2019-01-06) ==================================== %

% ----------------------------------------------------------------------- %
% -- Initializing mean measurements variables --------------------------- %
mean_error_mat_LMS_comp = [];
mean_MSD_mat_LMS_comp = [];
mean_orig_MSD_mat_LMS_comp = [];
update_counter_vec_LMS_comp = [];

% ----------------------------------------------------------------------- %
% -- Basic Simulation Settings ------------------------------------------ %
numberSamples = 12000;

mu_vec = [ 0.2 0.5 1.0 0.2 0.5 1.0 ];   % mu_L
alg_selected = 1;                       % LMS

%%  ==================================================================== %%
% ======================================================================= %
% == Running Simulation and Obtaining Results (II) ====================== %

for int_counter = 1:length(mu_vec) 
    int_counter
    
    % ------------------------------------------------------------------- %
    % ------------------------------------------------------------------- %
    [ mean_error_mat_LMS_comp(:,int_counter), mean_MSD_mat_LMS_comp(:,int_counter), ...
      mean_orig_MSD_mat_LMS_comp(:,int_counter), ~, update_counter_vec_LMS_comp(int_counter), ~ ] = ...
      graph_signal_reconstruction( graph_signal , bandlimited_graph , [ 1 1 ], ...
            D_s, U_f, diag(variance_vector_matrix(:,int_counter)), kappa_vec(int_counter), ...
            mu_vec(int_counter), alg_selected, DS_strategy_vec(int_counter), ...
            numberSamples, numbCheckedUpdates, ensemble, graph_signal_type) ;
    % ------------------------------------------------------------------- %
    % ------------------------------------------------------------------- %
end

%%  ==================================================================== %%
% ======================================================================= %
% == Saving important data (II) for plotting results in another script == %
save('./Simu_Results/02-simu-DS-LMS','mean_MSD_mat_LMS_comp','mean_error_mat_LMS_comp',...
    'mu_vec','variance_vector_matrix','ensemble');


%%  ==================================================================== %%
% ======================================================================= %
% == III - RLS Algorithm (2019-01-06) =================================== %

% ----------------------------------------------------------------------- %
% -- Initializing mean measurements variables --------------------------- %
mean_error_mat_RLS_comp = [];
mean_MSD_mat_RLS_comp = [];
mean_orig_MSD_mat_RLS_comp = [];
update_counter_vec_RLS_comp = [];

% ----------------------------------------------------------------------- %
% -- Basic Simulation Settings ------------------------------------------ %
numberSamples = 1500; 

beta_vec = [ 0.95 0.90 0.75 0.95 0.90 0.75 ];   % beta_R
alg_selected = 2;                               % RLS


%%  ==================================================================== %%
% ======================================================================= %
% == Running Simulation and Obtaining Results (III) ===================== %

for int_counter = 1:length(beta_vec) 
    int_counter
    
    % ------------------------------------------------------------------- %
    % ------------------------------------------------------------------- %
    [ mean_error_mat_RLS_comp(:,int_counter), mean_MSD_mat_RLS_comp(:,int_counter), ...
      mean_orig_MSD_mat_RLS_comp(:,int_counter), ~, update_counter_vec_RLS_comp(int_counter), ~ ] = ...
      graph_signal_reconstruction( graph_signal , bandlimited_graph , [ 1 1 ], ...
            D_s, U_f, diag(variance_vector_matrix(:,int_counter)), kappa_vec(int_counter), ...
            beta_vec(int_counter), alg_selected, DS_strategy_vec(int_counter), ...
            numberSamples, numbCheckedUpdates, ensemble, graph_signal_type) ;
    % ------------------------------------------------------------------- %
    % ------------------------------------------------------------------- %
end


%%  ==================================================================== %%
% ======================================================================= %
% == Saving important data (III) for plotting results in another script = %
save('./Simu_Results/02-simu-DS-RLS','mean_MSD_mat_RLS_comp','mean_error_mat_RLS_comp',...
    'beta_vec','variance_vector_matrix','ensemble');

% == END OF SCRIPT ====================================================== %
% ======================================================================= %
