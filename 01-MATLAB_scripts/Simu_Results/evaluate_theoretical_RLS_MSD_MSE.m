% ======================================================================= %
% ======================================================================= %

%clc;clear;

load('graph_data_1')

% ======================================================================= %
M = size(U_f,2) + 10; %10; %50;
% ----------------------------------------------------------------------- %
% -- Uncomment this part if it is necessary to find the indices to be
% -- sampled
% p = diag( eig_sampling_strategy(M, U_f) );
% save('p_selection','p')
% disp('p vector has been selected')
% ----------------------------------------------------------------------- %
load('p_selection')
D = diag(p);
% ======================================================================= %

%% ===================================================================== %%
% == Estimating Stationary Values for the RLS Algorithm ================= %

Cw_choice = 3;

switch Cw_choice
    case 1
        sigma_w = sqrt(0.01); %sqrt(0.01);
        variance_vector = sigma_w^2 *ones(299,1);
    case 2
        rng(1);
        variance_vector = 0.003 + 0.007*rand(299,1);
    case 3
        rng(2);
        variance_vector = 0.005 + 0.010*rand(299,1);
end

C_w = diag(variance_vector);

beta_vector = [ 0.95 0.75 ] ; %0.9 0.75 0.5 ];

%%

inv_mat = inv(U_f'*D*inv(C_w)*D*U_f);

S_r = inv_mat * ( U_f'*D*inv(C_w)*D*C_w*D*inv(C_w)*U_f ) * inv_mat;

for i = 1:length(beta_vector)

    beta = beta_vector(i);
    
    for counter = 1:299
        if (D(counter,counter) == 1)
            m = U_f(counter,:)';
            est_var_vec(counter,1) = C_w(counter,counter) + (1-beta)/(1+beta) * m' * S_r * m;
        else
            est_var_vec(counter,1) = 0;
        end
    end

    %figure;stem(est_var_vec)

    sum(est_var_vec)
                pause
               
    
    MSD_vec(i) = (1-beta)/(1+beta) * trace( S_r );
    MSE_vec(i) = sum(est_var_vec);
    
    %my_gamma_vec = sqrt(est_var_vec);
    
end

disp('=================================================================== ')
disp('MSD Values')
for int_counter = 1:length(beta_vector)
    disp( [ 'mu = ' num2str(beta_vector(int_counter)) ' - MSD = ' num2str( MSD_vec(int_counter) ) ...
        ' - MSD[dB] = ' num2str( 10*log10(MSD_vec(int_counter)) ) ] )
end
disp('=================================================================== ')
disp('=================================================================== ')
disp('MSE Values')
for int_counter = 1:length(beta_vector)
    disp( [ 'mu = ' num2str(beta_vector(int_counter)) ' - MSE = ' num2str( MSE_vec(int_counter) ) ...
        ' - MSE[dB] = ' num2str( 10*log10(MSE_vec(int_counter)) ) ] )
end
disp('=================================================================== ')
