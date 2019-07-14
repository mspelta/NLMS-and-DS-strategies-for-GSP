function [ threshold_vector ] = evaluate_threshold_vector(U_f, D, alg_param, algorithm_selection, C_w, DS_strategy)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

    N = size(U_f,1);

    switch algorithm_selection
        case 1  % --> LMS Algorithm
            % ----------------------------------------------------------- %
            
            mu = alg_param;
    
            A = U_f'*D*U_f;
            B = U_f'*D*C_w*D*U_f;

            S_n = zeros(size(A));

            for counter = 1:100   % I have used 100 for guaranteeing closeenough approximation
                C = mu*(B + A*S_n*A);
                S_n = sylvester(A,A,C);
            end

            for counter = 1:299
                if (D(counter,counter) == 1)
                    m = U_f(counter,:)';
                    var_vec(counter,1) = C_w(counter,counter) +  m' * S_n * m;
                else
                    var_vec(counter,1) = 0;
                end
            end
            % ----------------------------------------------------------- %
        case 2  % --> RLS Algorithm
            % ----------------------------------------------------------- %
            
            inv_mat = inv(U_f'*D*inv(C_w)*D*U_f);

            S_r = inv_mat * ( U_f'*D*inv(C_w)*D*C_w*D*inv(C_w)*U_f ) * inv_mat;

            beta = alg_param;
    
            for counter = 1:299
                if (D(counter,counter) == 1)
                    m = U_f(counter,:)';
                    var_vec(counter,1) = C_w(counter,counter) + (1-beta)/(1+beta) * m' * S_r * m;
                else
                    var_vec(counter,1) = 0;
                end
            end
 
            % ----------------------------------------------------------- %
            case 3 % --> NLMS Algorithm
            % ----------------------------------------------------------- %
            big_mat = inv(U_f'*D*U_f);
            err_sigma_vec = [];

            mu = alg_param;
            
            for i = 1:N
                u_i = U_f(i,:)';

                if ( D(i,i) == 1)
                    var_vec(i) = ( C_w(i,i)  + mu/(2-mu) * u_i' * big_mat * (U_f'*D*C_w*D*U_f) * big_mat * u_i ) ;
                else
                    var_vec(i) = 0;
                end

            end

            % ----------------------------------------------------------- %
    end

     if(DS_strategy == 1) % component-wise choice
        threshold_vector = sqrt(var_vec);
     else                 % l2-norm choice
        mean_var = sum(var_vec)/( sum(diag(D)) );
        threshold_vector = sum(var_vec)*ones(size(U_f,1),1);
     end
    
end

