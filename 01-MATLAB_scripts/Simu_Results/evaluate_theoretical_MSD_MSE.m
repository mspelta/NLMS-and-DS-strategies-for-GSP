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
var_v = 0.01;
mu_vector = [0.002 0.01 0.05 0.25 1.0];

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

    MSD_vec(i) = mu/(2-mu)*trace(inv(U_f'*D*U_f))
    MSE_vec(i) = sum(est_var_vec);
    
    %my_gamma_vec = sqrt(est_var_vec);
    
end

disp('=================================================================== ')
disp('MSD Values')
for int_counter = 1:length(MSD_vec)
    disp( [ 'mu = ' num2str(mu_vector(int_counter)) ' - MSD = ' num2str( 10*log10(MSD_vec(int_counter)) ) ] )
end
disp('=================================================================== ')
disp('=================================================================== ')
disp('MSE Values')
for int_counter = 1:length(MSE_vec)
    disp( [ 'mu = ' num2str(mu_vector(int_counter)) ' - MSE = ' num2str( 10*log10(MSE_vec(int_counter)) ) ] )
end
disp('=================================================================== ')
disp('=================================================================== ')
disp('Convergence Speed Values')
for int_counter = 1:length(MSE_vec)
    disp( [ 'mu = ' num2str(mu_vector(int_counter)) ' - Ang. Coeff. = ' num2str( log10(1 - mu_vector(int_counter)) ) ] )
end
disp('=================================================================== ')

%% =================================

%inv_mat_2 = inv(2*eye(200) - U_f'*D*U_f);
mu_vector = [0.5 0.6 0.7 0.8]; %[0.1 0.2 0.3 0.4]; %[0.02 0.03 0.04 0.05]; %[ 0.025 0.05 0.1 0.2 ]; % 0.25 0.5 ]; %0.25 0.5 1.0 2.0]; %[0.002 0.01 0.05 0.25 1.0 2.0];
A = U_f'*D*U_f; %single(U_f'*D*U_f);
var_v = 0.01;

MSD_vec_LMS = [];
MSE_vec_LMS = [];

for i = 1:length(mu_vector)

    mu = mu_vector(i);
    
    C = mu*var_v*U_f'*D*U_f;
    inv_mat_2 = sylvester(A,A,C);
    
    %pause
    
    for counter = 1:299
        if (D(counter,counter) == 1)
            m = U_f(counter,:)';
            est_var_vec_LMS(counter,1) = var_v +  m' * inv_mat_2 * m;
        else
            est_var_vec_LMS(counter,1) = 0;
        end
    end

    %figure;stem(est_var_vec_LMS)
    %pause
    
    MSD_vec_LMS(i) = trace( (eye(200) - mu*A)*inv_mat_2*(eye(200) - mu*A) ) + mu^2 * var_v * trace( A ) ;
    MSE_vec_LMS(i) = sum(est_var_vec_LMS);
    
    %my_gamma_vec = sqrt(est_var_vec);
    
end

disp('=================================================================== ')
disp('MSE Values')
for int_counter = 1:length(MSE_vec_LMS)
    disp( [ 'mu = ' num2str(mu_vector(int_counter)) ' - MSD = ' num2str( MSD_vec_LMS(int_counter) ) ' - MSD[dB] = ' num2str( 10*log10(MSD_vec_LMS(int_counter)) ) ] )
end
disp('=================================================================== ')
disp('MSE Values')
for int_counter = 1:length(MSE_vec_LMS)
    disp( [ 'mu = ' num2str(mu_vector(int_counter)) ' - MSE = ' num2str( MSE_vec_LMS(int_counter) ) ' - MSE[dB] = ' num2str( 10*log10(MSE_vec_LMS(int_counter)) ) ] )
end
disp('=================================================================== ')

%% =================================
% Diagonal C_w with different variance values (General case!!)

%inv_mat_2 = inv(2*eye(200) - U_f'*D*U_f);
mu_vector = [ 0.05 ]; %[ 0.02 0.03 0.04 0.05 ]; %[0.02 0.05 0.10 0.15 0.20 0.25]; %[0.5 0.6 0.7 0.8]; %[0.1 0.2 0.3 0.4]; %[0.02 0.03 0.04 0.05]; %[ 0.025 0.05 0.1 0.2 ]; % 0.25 0.5 ]; %0.25 0.5 1.0 2.0]; %[0.002 0.01 0.05 0.25 1.0 2.0];
A = U_f'*D*U_f; %single(U_f'*D*U_f);
%var_v = 0.01;

rng(1);
variance_vector = 0.003 + 0.007*rand(299,1);
C_w = diag(variance_vector);

MSD_vec_LMS = [];
MSE_vec_LMS = [];



for i = 1:length(mu_vector)

    mu = mu_vector(i);
    
    C = mu*U_f'*D*C_w*D*U_f;
    inv_mat_2 = sylvester(A,A,C);
    
    %pause
    
    for counter = 1:299
        if (D(counter,counter) == 1)
            m = U_f(counter,:)';
            est_var_vec_LMS(counter,1) = variance_vector(counter) +  m' * inv_mat_2 * m;   % Equation (23) and (24)
        else
            est_var_vec_LMS(counter,1) = 0;
        end
    end

    %figure;stem(est_var_vec_LMS)
    %pause
    
    MSD_vec_LMS(i) = trace( (eye(200) - mu*A)*inv_mat_2*(eye(200) - mu*A) ) + mu^2 * trace( C ) ;
    MSE_vec_LMS(i) = sum(est_var_vec_LMS);
    
    %my_gamma_vec = sqrt(est_var_vec);
    
end

disp('=================================================================== ')
disp('MSE Values')
for int_counter = 1:length(MSE_vec_LMS)
    disp( [ 'mu = ' num2str(mu_vector(int_counter)) ' - MSD = ' num2str( MSD_vec_LMS(int_counter) ) ' - MSD[dB] = ' num2str( 10*log10(MSD_vec_LMS(int_counter)) ) ] )
end
disp('=================================================================== ')
disp('MSE Values')
for int_counter = 1:length(MSE_vec_LMS)
    disp( [ 'mu = ' num2str(mu_vector(int_counter)) ' - MSE = ' num2str( MSE_vec_LMS(int_counter) ) ' - MSE[dB] = ' num2str( 10*log10(MSE_vec_LMS(int_counter)) ) ] )
end
disp('=================================================================== ')

%% =================================
% Diagonal C_w with different variance values (General case!!)

%inv_mat_2 = inv(2*eye(200) - U_f'*D*U_f);
mu_vector = [ 0.05 0.075 0.1]; %[ 0.02 0.03 0.04 0.05 ]; %[0.02 0.05 0.10 0.15 0.20 0.25]; %[0.5 0.6 0.7 0.8]; %[0.1 0.2 0.3 0.4]; %[0.02 0.03 0.04 0.05]; %[ 0.025 0.05 0.1 0.2 ]; % 0.25 0.5 ]; %0.25 0.5 1.0 2.0]; %[0.002 0.01 0.05 0.25 1.0 2.0];
A = U_f'*D*U_f; %single(U_f'*D*U_f);
%var_v = 0.01;

rng(1);
variance_vector = 0.003 + 0.007*rand(299,1);
%variance_vector = 0.01*ones(299,1);
C_w = diag(variance_vector);

MSD_vec_LMS = [];
MSE_vec_LMS = [];



for i = 1:length(mu_vector)

    mu = mu_vector(i);
    
    C = mu*U_f'*D*C_w*D*U_f;
    inv_mat_2 = sylvester(A,A,C);
    
    %pause
    
    for counter = 1:299
        if (D(counter,counter) == 1)
            m = U_f(counter,:)';
            est_var_vec_LMS(counter,1) = variance_vector(counter) +  m' * inv_mat_2 * m;   % Equation (23) and (24)
        else
            est_var_vec_LMS(counter,1) = 0;
        end
    end

    %figure;stem(est_var_vec_LMS)
    %pause
    
    MSD_vec_LMS(i) = trace( (eye(200) - mu*A)*inv_mat_2*(eye(200) - mu*A) ) + mu * trace( C ) ;
    MSE_vec_LMS(i) = sum(est_var_vec_LMS);
    
    %my_gamma_vec = sqrt(est_var_vec);
    
end

disp('=================================================================== ')
disp('MSD Values')
for int_counter = 1:length(MSE_vec_LMS)
    disp( [ 'mu = ' num2str(mu_vector(int_counter)) ' - MSD = ' num2str( MSD_vec_LMS(int_counter) ) ' - MSD[dB] = ' num2str( 10*log10(MSD_vec_LMS(int_counter)) ) ] )
end
disp('=================================================================== ')
disp('MSE Values')
for int_counter = 1:length(MSE_vec_LMS)
    disp( [ 'mu = ' num2str(mu_vector(int_counter)) ' - MSE = ' num2str( MSE_vec_LMS(int_counter) ) ' - MSE[dB] = ' num2str( 10*log10(MSE_vec_LMS(int_counter)) ) ] )
end
disp('=================================================================== ')

%%

prev_S = zeros(200);

for counter = 1:5
    S = mu*A*prev_S*A + mu*C;

    errorMat = prev_S*A + A*prev_S - C;

    errorMat(1:10,1:10)
    
    prev_S = S;
    pause
end