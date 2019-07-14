% ======================================================================= %
% ======================================================================= %
% == Script: graph_signal_reconstruction
% == Responsible: Marcelo Jorge Mendes Spelta - Date: 2018/12/15
% == E-mail: marcelojorgespelta@poli.ufrj.br / marcelo.spelta@smt.ufrj.br
% ======================================================================= %
% == Short Description: This auxiliar general script is used for recons-
% == tructing the original bandlimited graph signal from a reduced number of
% == noisy measurements using the DS-LMS, DS-RLS and DS-NLMS adaptive 
% == reconstruction strategies.
% ======================================================================= %
% ======================================================================= %
    
function [ mean_MSE_vector, mean_MSD_vector, mean_x_vec_matrix, ...
           mean_update_counter, mean_elapsedTime_vector ] = ...
            graph_signal_reconstruction( x_ref , ...
            D, U_f, C_w, kappa, alg_factor, alg_selection, DS_strategy, numberSamples, ...
            numbCheckedUpdates, ensemble )

% if graph_signal_type == 1
%     x_signal = bandlimited_graph; % reference bandlimited graph signal x_o
% else
%     x_signal = original_x_signal;
% end
% 
% 
% 
% % ======================================================================= %
% % == Generating x_o[k] and its scalar variations for all the iterations = %
% % == considered.
% % --------------------------------------
% for counter = 1:(numberSamples/2)
%     x_ref(:,counter) = amplitude_values(1) * x_signal;
% end
% for counter = (numberSamples/2 + 1):numberSamples
%     x_ref(:,counter) = amplitude_values(2) * x_signal; %1.1*x_signal;
% end
% --------------------------------------

N = size(x_ref,1);   % total number of vertices for the used graph signal
F = size(U_f, 2);       % amount of frequency components considered
%D = diag(p);            % Sampling matrix D_s

% ======================================================================= %
% ======================================================================= %

% ----------------------------------------------------------------------- %
% -- Initialization of Figure of Merit (FoM) auxiliar matrices and other
% -- auxiliar variables usid in this script
MSE_vector_matrix = []; MSD_vector_matrix = []; elapsedTime_vector_matrix = [];
x_vec_matrix = [];
update_counter_vec = [];
%C_w = diag( variance_vector );  % Covariance matrix
delta = 1e-3;                   % small parameter used in the RLS algorithm for initializing the PSI matrix

% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% -- Obtaining the reference constraint parameters for the data-selection 
% -- strategies. If (DS_strategy == 1), the ref_gamma_vec is a vector
% -- evaluated according to the component-wise error constraint formulas.
% -- On the other hand, if (DS_strategy == 2), the constraint parameter is
% -- a scalar value, which is taken as the first component of the vector
% -- ref_gamma_vec, i.e., ref_gamma_vec(1). 
%ref_gamma_vec = kappa * evaluate_threshold_vector(U_f, D, alg_factor, alg_selection, C_w, DS_strategy);
var_vec = evaluate_variance_vector( U_f, D, alg_factor, alg_selection, C_w );
% ----------------------------------------------------------------------- %
% ======================================================================= %
% ======================================================================= %
% -- Ensemble Loop ------------------------------------------------------ %

% --------------
switch(alg_selection)
    case 1      % LMS algorithm
        B_l = U_f * U_f';                           
        mu_Bl_matrix = alg_factor * B_l;
    case 2      % RLS algorithm
        B_r = U_f' * D * inv(C_w);
        M_prime = U_f' * D * inv(C_w) * D * U_f;  
    case 3      % NLMS algorithm
        B_n = U_f * inv( U_f' * D * U_f ) * U_f'; 
        mu_Bn_matrix = alg_factor * B_n;
end
% --------------------

x_vec_tensor = [];

for extCounter = 1:ensemble
    % ------------------------------------------------------------------- %
    % -- Initialization of internal auxiliar variables
    MSE_vector = []; MSD_vector = []; elapsedTime_vector = [];
    x_vec = zeros(N,1);
    PSI = delta* eye( F ); %psi = zeros( F, 1); % RLS algorithm internal variables
    update_counter = 0;
    % ------------------------------------------------------------------- %
    x_vec_matrix = [];
    for itCounter = 1:numberSamples
        
%         max(abs(x_ref(:,itCounter) - bandlimited_graph))
        %a_aux = x_ref(:,itCounter) - bandlimited_graph;
%         disp('Key')
%         pause
        % --------------------------------------------------------------- %
        % -- Generating sampled noisy measurements ---------------------- %
        x_w = D* ( x_ref(:,itCounter) + sqrt( diag(C_w) ).*randn(N,1) );
        % --------------------------------------------------------------- %
        
        % --------------------------------------------------------------- %
        % -- Error Signal ----------------------------------------------- %
        error_signal = D*(x_w - x_vec); % Current error signal vector e[k]
%         a_aux(295:299)
%         error_signal(295:299)
%         pause
        
        MSE_vector = [MSE_vector;
                        norm( error_signal, 2 )^2 ]; % for plotting later!
        % --------------------------------------------------------------- %
        
        % --------------------------------------------------------------- %
        % -- Mean-Squared Deviation Metrics ----------------------------- %
        MSD = norm(x_ref(:,itCounter) - x_vec, 2)^2;   
        MSD_vector = [MSD_vector;
                      MSD];
        
        % --------------------------------------------------------------- %
        
        % =============================================================== %
        % =============================================================== %
        % == Implementation of the different adaptive algorithms updates  %
        switch(alg_selection)
            % ----------------------------------------------------------- %
            % -- Original LMS ------------------------------------------- %
            case 1
                [ x_vec, update_flag, elapsedTime ] = ...
                    DS_LMS_GSP( x_w, x_vec, D, mu_Bl_matrix, var_vec, DS_strategy, kappa );
            % ----------------------------------------------------------- %
            % -- Original RLS ------------------------------------------- %
            case 2
                [ x_vec, PSI, update_flag, elapsedTime ] = ...
                    DS_RLS_GSP( x_w, x_vec, U_f, D, M_prime, B_r, alg_factor, PSI, var_vec, DS_strategy, kappa );
            % ----------------------------------------------------------- %
            % ----------------------------------------------------------- %
            % -- DS-NLMS Algorithm (20/11/2018) ------------------------ %
            case 3
                [ x_vec, update_flag, elapsedTime ] = ...
                    DS_NLMS_GSP( x_w, x_vec, D, mu_Bn_matrix, var_vec, DS_strategy, kappa );
            % ----------------------------------------------------------- %
            % ----------------------------------------------------------- %
        end
        % =============================================================== %
        % =============================================================== %
        x_vec_matrix = [x_vec_matrix x_vec];
        % --------------------------------------------------------------- %
        elapsedTime_vector = [ elapsedTime_vector ;
                               elapsedTime ];
        % --------------------------------------------------------------- %
        % -- Checks if the update_counter variable must be incremented or not
        if ( (update_flag == true) && (itCounter > (numberSamples - numbCheckedUpdates) ) )
            update_counter = update_counter + 1;
        end
        % --------------------------------------------------------------- %
    end
    
    % ------------------------------------------------------------------- %
    % -- Storing important information about each algorithm run in order 
    % -- to perform an average operation later. 
    MSE_vector_matrix = [MSE_vector_matrix MSE_vector];
    MSD_vector_matrix = [MSD_vector_matrix MSD_vector];
    x_vec_tensor(:,:,extCounter) = x_vec_matrix;
    elapsedTime_vector_matrix = [ elapsedTime_vector_matrix elapsedTime_vector ];
    update_counter_vec = [ update_counter_vec update_counter/numbCheckedUpdates ];
    % ------------------------------------------------------------------- %
%     size(x_ref)
%     plot(1:95,x_vec_matrix(100,:),'b',1:95,x_ref(100,:),'r')
%     pause
        
    extCounter   % indicates to the user the last algorithm run performed
    
end

% ======================================================================= %
% ----------------------------------------------------------------------- %
% -- Computing ensemble average values that must be returned to the user %
mean_MSE_vector = [];
mean_MSD_vector = [];
mean_elapsedTime_vector = [];
mean_update_counter = mean(update_counter_vec);
mean_x_vec_matrix = [];

for counter = 1:numberSamples
    mean_MSE_vector(counter) = mean( MSE_vector_matrix(counter,:) )';
    mean_MSD_vector(counter) = mean( MSD_vector_matrix(counter,:) )';
    mean_elapsedTime_vector(counter) = mean( elapsedTime_vector_matrix(counter,:) )';
%     if (counter <= length(x_vec) )
%         mean_x_vec(counter) = mean( x_vec_matrix(counter,:) )';
%     end
end
for i = 1:205
    for j = 1:95
        mean_x_vec_matrix(i,j) = mean(x_vec_tensor(i,j,:));
    end
end
% plot(1:95,mean_x_vec_matrix(100,:),'b--',1:95,x_ref(100,:),'r--')
% pause
% ======================================================================= %
% ======================================================================= %

end
% == END OF FUNCTION ==================================================== %
% ======================================================================= %