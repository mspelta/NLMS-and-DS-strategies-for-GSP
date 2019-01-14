% ======================================================================= %
% ======================================================================= %
% == Script: script03_DataSelection_final =============================== %
% == Responsible: Marcelo Jorge Mendes Spelta - Date: 2019/01/12 ======== %
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

%% ===================================================================== %%
% == Initializing general variables used in all simulations ============= %

ensemble = 1000; %1000;  
numbCheckedUpdates = 1000;

%%  ==================================================================== %%
% ======================================================================= %
% == I - DS-NLMS Algorithm ============================================== %

% ----------------------------------------------------------------------- %
% -- Initializing mean measurements variables --------------------------- %
mean_error_mat_LMS_comp = [];
mean_MSD_mat_LMS_comp = [];
mean_orig_MSD_mat_LMS_comp = [];
update_counter_vec_LMS_comp = [];

% ----------------------------------------------------------------------- %
% -- Basic Simulation Settings ------------------------------------------ %
numberSamples = 2500;

mu = 0.1;                   % mu_N 
alg_selection = 3;          % NLMS 
DS_strategy_vec = [ 1*ones(1,6) 2*ones(1,6) 1*ones(1,6) 2*ones(1,6) ]; 
kappa_vec = [ 2.50 2.75 3.00 3.25 3.50 3.75 0.90 0.95 1.00 1.05 1.10 1.15 ... 
              2.50 2.75 3.00 3.25 3.50 3.75 0.90 0.95 1.00 1.05 1.10 1.15 ];

% ----------------------------------------------------------------------- %
% -- Generating diagonal variance vectors. That will be used to form C_w  %
variance_vector_matrix = 0.01*ones(299,12);
rng(2);
variance_vector = 0.005 + 0.010*rand(299,1);   % average is about variance 0.01 
variance_vector_matrix = [ variance_vector_matrix variance_vector*ones(1,12) ];
% ----------------------------------------------------------------------- %

% ======================================================================= %  

%%  ==================================================================== %%
% ======================================================================= %
% == Running Simulation and Obtaining Results (I) ======================= %

for int_counter = 1:length(kappa_vec) 
    int_counter
    % =================================================================== %
    % =================================================================== %
    [ mean_error_mat_LMS_comp(:,int_counter), mean_MSD_mat_LMS_comp(:,int_counter), ...
      mean_orig_MSD_mat_LMS_comp(:,int_counter), ~, ...
      update_counter_vec_LMS_comp(int_counter), ~ ] = ...
      graph_signal_reconstruction( graph_signal , bandlimited_graph , [ 1 1 ], ...
            D, U_f, diag(variance_vector_matrix(:,int_counter)), kappa_vec(int_counter), ...
            mu, alg_selection, DS_strategy_vec(int_counter), numberSamples, numbCheckedUpdates, ensemble) ;
    % =================================================================== %
    % =================================================================== %
end

%%  ==================================================================== %%
% ======================================================================= %
% == Saving important data (I) for plotting results in another script === %
save('./Simu_Results/03-simu-DS-NLMS_DataSel','mean_MSD_mat_LMS_comp','update_counter_vec_LMS_comp',...
    'mu','variance_vector_matrix','DS_strategy_vec','kappa_vec','ensemble','numbCheckedUpdates');


%%  ==================================================================== %%
% ======================================================================= %
% == II - DS-LMS Algorithm ============================================== %

% ----------------------------------------------------------------------- %
% -- Initializing mean measurements variables --------------------------- %
mean_error_mat_LMS_comp = [];
mean_MSD_mat_LMS_comp = [];
mean_orig_MSD_mat_LMS_comp = [];
update_counter_vec_LMS_comp = [];

% ----------------------------------------------------------------------- %
% -- Basic Simulation Settings ------------------------------------------ %
numberSamples = 10000;

mu = 0.5;                   % mu_L
alg_selection = 1;          % LMS 
DS_strategy_vec = [ 1*ones(1,4) 2*ones(1,4) 1*ones(1,4) 2*ones(1,4) ]; 
kappa_vec = [ 2.50 3.00 3.50 3.75 0.90 1.00 1.10 1.15 ... 
              2.50 3.00 3.50 3.75 0.90 1.00 1.10 1.15 ];
          
% ----------------------------------------------------------------------- %
% -- Generating diagonal variance vectors. That will be used to form C_w  %       
variance_vector_matrix = [ 0.01*ones(299,8) ];
rng(2);
variance_vector = 0.005 + 0.010*rand(299,1);
variance_vector_matrix = [ variance_vector_matrix  variance_vector*ones(1,8) ];
% ----------------------------------------------------------------------- %

%%  ==================================================================== %%
% ======================================================================= %
% == Running Simulation and Obtaining Results (II) ====================== %

for int_counter = 1:length(kappa_vec)  
    int_counter
    % =================================================================== %
    % =================================================================== %
    [ mean_error_mat_LMS_comp(:,int_counter), mean_MSD_mat_LMS_comp(:,int_counter), ...
      mean_orig_MSD_mat_LMS_comp(:,int_counter), ~, ...
      update_counter_vec_LMS_comp(int_counter), ~ ] = ...
      graph_signal_reconstruction( graph_signal , bandlimited_graph , [ 1 1 ], ...
            D, U_f, diag(variance_vector_matrix(:,int_counter)), kappa_vec(int_counter), ...
            mu, alg_selection, DS_strategy_vec(int_counter), numberSamples, numbCheckedUpdates, ensemble) ;
    % =================================================================== %
    % =================================================================== %
end


%%  ==================================================================== %%
% ======================================================================= %
% == Saving important data (II) for plotting results in another script == %

save('./Simu_Results/03-simu-DS-LMS','mean_MSD_mat_LMS_comp','update_counter_vec_LMS_comp',...
    'mu','variance_vector_matrix','DS_strategy_vec','kappa_vec','ensemble','numbCheckedUpdates');


%%  ==================================================================== %%
% ======================================================================= %
% == III - DS-RLS Algorithm ============================================= %

% ----------------------------------------------------------------------- %
% -- Initializing mean measurements variables --------------------------- %
mean_error_mat_RLS_comp = [];
mean_MSD_mat_RLS_comp = [];
mean_orig_MSD_mat_RLS_comp = [];
update_counter_vec_RLS_comp = [];

% ----------------------------------------------------------------------- %
% -- Basic Simulation Settings ------------------------------------------ %

numberSamples = 2000; 

beta = 0.9;                         % beta_R
alg_selection = 2;                  % RLS 
DS_strategy_vec = [ 1*ones(1,4) 2*ones(1,4) 1*ones(1,4) 2*ones(1,4) ]; 
kappa_vec = [ 2.50 3.00 3.50 3.75 0.90 1.00 1.10 1.15 ... 
              2.50 3.00 3.50 3.75 0.90 1.00 1.10 1.15 ];

% ----------------------------------------------------------------------- %
% -- Generating diagonal variance vectors. That will be used to form C_w  % 
variance_vector_matrix = [ 0.01*ones(299,8) ];
rng(2);
variance_vector = 0.005 + 0.010*rand(299,1);
variance_vector_matrix = [ variance_vector_matrix variance_vector*ones(1,8) ];
% ----------------------------------------------------------------------- %


%%  ==================================================================== %%
% ======================================================================= %
% == Running Simulation and Obtaining Results (III) ===================== %

for int_counter = 1:length(kappa_vec) 
    int_counter
    % =================================================================== %
    % =================================================================== %
    [ mean_error_mat_RLS_comp(:,int_counter), mean_MSD_mat_RLS_comp(:,int_counter), ...
      mean_orig_MSD_mat_RLS_comp(:,int_counter), ~, ...
      update_counter_vec_RLS_comp(int_counter), ~ ] = ...
      graph_signal_reconstruction( graph_signal , bandlimited_graph , [ 1 1 ], ...
            D, U_f, diag(variance_vector_matrix(:,int_counter)), kappa_vec(int_counter), ...
            beta, alg_selection, DS_strategy_vec(int_counter), numberSamples, numbCheckedUpdates, ensemble) ;
    % =================================================================== %
    % =================================================================== %
end

%%  ==================================================================== %%
% ======================================================================= %
% == Saving important data (III) for plotting results in another script = %

save('./Simu_Results/03-simu-DS-RLS','mean_MSD_mat_RLS_comp','update_counter_vec_RLS_comp',...
    'beta','variance_vector_matrix','DS_strategy_vec','kappa_vec','ensemble','numbCheckedUpdates');

% ======================================================================= %
% ======================================================================= %
