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


%% ===================================================================== %%
% == Estimating Stationary Values for the RLS Algorithm ================= %

mu_vector = [];
mu_vector = [ 0.1 0.2 0.5 ];

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

%inv_mat = inv(U_f'*D*U_f);

%%

for i = 1:length(mu_vector)

    mu = mu_vector(i);
    
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

    MSD_vec(i) = trace(S_n);
    MSE_vec(i) = sum(est_var_vec);
    
    %my_gamma_vec = sqrt(est_var_vec);
    
end

disp('================================')
disp('Theoretical Results')
disp('=================================================================== ')
disp('MSD Values')
for int_counter = 1:length(mu_vector)
    disp( [ 'mu = ' num2str(mu_vector(int_counter)) ' - MSD = ' num2str( MSD_vec(int_counter) ) ...
        ' - MSD[dB] = ' num2str( 10*log10(MSD_vec(int_counter)) ) ] )
end
disp('=================================================================== ')
disp('=================================================================== ')
disp('MSE Values')
for int_counter = 1:length(mu_vector)
    disp( [ 'mu = ' num2str(mu_vector(int_counter)) ' - MSE = ' num2str( MSE_vec(int_counter) ) ...
        ' - MSE[dB] = ' num2str( 10*log10(MSE_vec(int_counter)) ) ] )
end
disp('=================================================================== ')
