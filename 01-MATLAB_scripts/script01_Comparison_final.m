% ======================================================================= %
% ======================================================================= %
% == Script: script01_Comparison_final ================================== %
% == Responsible: Marcelo Jorge Mendes Spelta - Date: 2019/01/06 ======== %
% == E-mail: marcelo.spelta@smt.ufrj.br ================================= %
% ======================================================================= %
% ======================================================================= %

clc;clear;

% ======================================================================= %
% -- Loading original graph signal, bandlimited graph signal, U_f and  -- %
% -- index of vertices to be sampled. Graph signal is obtained from ----- %
% -- the INMET dataset. ------------------------------------------------- %
load('scenario_data')
D = diag(S);        % D_s -> Sampling matrix
% ======================================================================= %

%%  ==================================================================== %%
% ======================================================================= %
% == I - Comparison among the LMS, RLS and NLMS Algorithm (2019-01-06) == %

% ----------------------------------------------------------------------- %
% -- Initializing mean measurements variables --------------------------- %
mean_error_mat_comp = [];
mean_MSD_mat_comp = [];
mean_elapsedTime_mat_comp = [];

% ----------------------------------------------------------------------- %
% -- Basic Simulation Settings ------------------------------------------ %
numberSamples = 5000; 
ensemble = 1000; 
numbCheckedUpdates = 1000;

% alg_param_vec = [ mu_L beta_R mu_N / mu_L beta_R mu_N / mu_L beta_R mu_N ]
alg_param_vec = [ 0.280 0.930 0.070 0.280 0.930 0.070 0.721 0.792 0.208 ]; 
% alg_selection_vec = [ LMS RLS NLMS / LMS RLS NLMS / LMS RLS NLMS ]
alg_selection_vec = [ 1 2 3 1 2 3 1 2 3 ];  
DS_strategy_vec = ones(1,9); % Not important since we use kappa=0
kappa_vec = zeros(1,9);      % kappa=0 so we don't use the DS strategies
% ----------------------------------------------------------------------- %
% -- Generating diagonal variance vectors. That will be used to form C_w % 
variance_vector_matrix = [ 0.001*ones(299,3) ]; 
rng(2);
variance_vector = 0.005 + 0.010*rand(299,1);
variance_vector_matrix = [ variance_vector_matrix variance_vector*ones(1,6) ];
% ----------------------------------------------------------------------- %
% ======================================================================= %

%%  ==================================================================== %%
% ======================================================================= %
% == Running Simulation and Obtaining Results =========================== %

for int_counter = 1:length(alg_param_vec) 
    int_counter
    % ------------------------------------------------------------------- %
    % ------------------------------------------------------------------- %
    [ mean_error_mat_comp(:,int_counter), mean_MSD_mat_comp(:,int_counter), ...
        ~, ~, ~, mean_elapsedTime_mat_comp(:,int_counter) ] = ...
      graph_signal_reconstruction( graph_signal , bandlimited_graph , [ 1 1.2 ], ... 
            D, U_f, diag(variance_vector_matrix(:,int_counter)), kappa_vec(int_counter), ...
            alg_param_vec(int_counter), alg_selection_vec(int_counter), DS_strategy_vec(int_counter), ...
            numberSamples, numbCheckedUpdates, ensemble) ;
    % ------------------------------------------------------------------- %
    % ------------------------------------------------------------------- %
end

% ======================================================================= %

%%  ==================================================================== %%
% ======================================================================= %
% == Saving important data for plotting results in another script ======= %

save('./Simu_Results/01-simu-Comparison_LMS-RLS-NLMS','mean_MSD_mat_comp',...
    'mean_elapsedTime_mat_comp','alg_param_vec','alg_selection_vec','ensemble');

% ======================================================================= %
% ======================================================================= %