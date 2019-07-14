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

inv_mat = inv(U_f'*D*U_f);
var_v = 0.1;
mu_vector = [ 0.05 0.1 0.25 0.5 ];

%%

for i = 1:length(mu_vector)

    mu = mu_vector(i);
    
    for counter = 1:299
        if (D(counter,counter) == 1)
            m = U_f(counter,:)';
            est_var_vec(counter,1) = var_v + var_v * mu^2 / ( 1 - (1 - mu)^2) * m' * inv_mat * m;
        else
            est_var_vec(counter,1) = 0;
        end
    end

    %figure;stem(est_var_vec)

    MSD_vec(i) = mu/(2-mu)*var_v*trace(inv(U_f'*D*U_f))
    MSE_vec(i) = sum(est_var_vec);
    
    %my_gamma_vec = sqrt(est_var_vec);
    
end

disp('=================================================================== ')
disp('MSD Values')
for int_counter = 1:length(MSD_vec)
    disp( [ 'mu = ' num2str(mu_vector(int_counter)) ' - MSD = ' num2str( MSD_vec(int_counter) ) ...
        ' - MSD[dB] = ' num2str( 10*log10(MSD_vec(int_counter)) ) ] )
end
disp('=================================================================== ')
disp('=================================================================== ')
disp('MSE Values')
for int_counter = 1:length(MSE_vec)
    disp( [ 'mu = ' num2str(mu_vector(int_counter)) ' - MSE = ' num2str( MSE_vec(int_counter) ) ...
        ' - MSE[dB] = ' num2str( 10*log10(MSE_vec(int_counter)) ) ] )
end
disp('=================================================================== ')
disp('=================================================================== ')
disp('Convergence Speed Values')
for int_counter = 1:length(MSE_vec)
    disp( [ 'mu = ' num2str(mu_vector(int_counter)) ' - Ang. Coeff. = ' num2str( log10(1 - mu_vector(int_counter)) ) ] )
end
disp('=================================================================== ')

%% ===================================================================== %%

inv_mat = inv(U_f'*D*U_f);
%var_v = 0.1;
rng(2);
variance_vector = 0.003 + 0.007*rand(299,1);
C_w = diag(variance_vector);
mu_vector = [ 0.05 0.1 0.25 0.5 ];

%%

for i = 1:length(mu_vector)

    mu = mu_vector(i);
    
    for counter = 1:299
        if (D(counter,counter) == 1)
            m = U_f(counter,:)';
            est_var_vec(counter,1) = C_w(counter,counter) + mu / ( 2 - mu) * m' * inv_mat * (U_f'*D*C_w*D*U_f) * inv_mat*  m;
        else
            est_var_vec(counter,1) = 0;
        end
    end

    %figure;stem(est_var_vec)

    MSD_vec(i) = mu/(2-mu)* trace( inv_mat * (U_f'*D*C_w*D*U_f) * inv_mat )
    MSE_vec(i) = sum(est_var_vec);
    
    %my_gamma_vec = sqrt(est_var_vec);
    
end

disp('=================================================================== ')
disp('MSD Values')
for int_counter = 1:length(MSD_vec)
    disp( [ 'mu = ' num2str(mu_vector(int_counter)) ' - MSD = ' num2str( MSD_vec(int_counter) ) ...
        ' - MSD[dB] = ' num2str( 10*log10(MSD_vec(int_counter)) ) ] )
end
disp('=================================================================== ')
disp('=================================================================== ')
disp('MSE Values')
for int_counter = 1:length(MSE_vec)
    disp( [ 'mu = ' num2str(mu_vector(int_counter)) ' - MSE = ' num2str( MSE_vec(int_counter) ) ...
        ' - MSE[dB] = ' num2str( 10*log10(MSE_vec(int_counter)) ) ] )
end
disp('=================================================================== ')
disp('=================================================================== ')
disp('Convergence Speed Values')
for int_counter = 1:length(MSE_vec)
    disp( [ 'mu = ' num2str(mu_vector(int_counter)) ' - Ang. Coeff. = ' num2str( log10(1 - mu_vector(int_counter)) ) ] )
end
disp('=================================================================== ')
