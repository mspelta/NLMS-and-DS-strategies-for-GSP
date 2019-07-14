% ======================================================================= %
% == COPPE/UFRJ - Programa de Engenharia Eletrica (PEE) ================= %
% == Function: eig_sampling_strategy ==================================== %
% == Responsible: Marcelo Jorge Mendes Spelta - Date: 2019/03/29 ======== %
% == E-mail: marcelo.spelta@smt.ufrj.br ================================= %
% ======================================================================= %

function [ D, positionsVector ] = eig_sampling_strategy( M, U_f )

    threshold = 1e-6;   % zero threshold used

    [numb_lines, ~] = size(U_f);
    total_indices = numb_lines;
    
    indices_vector = zeros(total_indices, 1);
    numb_chosen_samples = 0;
    
    while ( numb_chosen_samples < M)
        min_nonzero_eig_vec = zeros(total_indices, 1);
        
        for counter = 1:total_indices
            %counter
            if (indices_vector(counter) == 0)
                indices_vector(counter) = 1;
                sorted_eigenvalue_list = sort( eig( U_f' * diag(indices_vector) * U_f ) );
                indices_vector(counter) = 0;
                
                exitLoop = false;
                intCounter = 1;
                while( (exitLoop == false) && (intCounter <= length( sorted_eigenvalue_list ) ) )
                    if ( sorted_eigenvalue_list(intCounter) > threshold )
                        exitLoop = true;
                        min_nonzero_eig_vec(counter) = sorted_eigenvalue_list(intCounter);
                    else
                        intCounter = intCounter + 1;
                    end
                end
            else
                min_nonzero_eig_vec(counter) = 0; % it is not taken into account
            end
            
        end
        
        [~, max_index] = max( min_nonzero_eig_vec );
        indices_vector(max_index) = 1;
        
         numb_chosen_samples =  numb_chosen_samples + 1
    end
    
    D = diag(indices_vector);
    
    positionsVector = [];
    for counter = 1:numb_lines
        if ( D(counter,counter) == 1 )
            positionsVector = [positionsVector counter];
        end
    end
    
end

% == END OF SCRIPT ====================================================== %
% ======================================================================= %