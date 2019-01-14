%% ===================================================================== %%
% ======================================================================= %
% == I - GSP NLMS Algorithm ============================================= %
% ======================================================================= %

clc; clear;

numbCheckedUpdates = 1000;

% ======================================================================= %
% -- Loading original graph signal, bandlimited graph signal, U_f and  -- %
% -- index of vertices to be sampled. Graph signal is obtained from ----- %
% -- the INMET dataset. ------------------------------------------------- %
load('../scenario_data')
D = diag(S);        % D_s -> Sampling matrix
% ======================================================================= %

% ----------------------------------------------
articleTable = zeros(10,13);
factors_vec = [ 0.05 0.10 0.25 0.50 0.2 0.5 1.0 0.95 0.9 0.75 ];
articleTable(:,1) = factors_vec;
% ----------------------------------------------

%% ===================================================================== %%

load('02-simu-DS-NLMS')

disp('=================================================================== ')
disp(['Ensemble of ' num2str(ensemble) ' runs'])
disp('Steady-state values taken considering the average of the last 2000 time instants.')
disp('=================================================================== ')
disp(['mu value = ' num2str(mu_vec(1))])
disp('=================================================================== ')
%% --------------------------------------------------------------------- %%
t = 1:length( mean_MSD_mat_LMS_comp(:,1) );
% ----------------------------------------------------------------------- %
% -- MSD (Mean Square Deviation) ---------------------------------------- %
figure
for int_counter = 1:length(mu_vec) 
    plot( t, 10*log10( mean_MSD_mat_LMS_comp(:,int_counter) ), 'LineWidth', 2 )
    hold on
end
ylabel('MSD [dB]','Interpreter','latex','fontsize',20)
xlabel('Number of Iterations $k$','Interpreter','latex','fontsize',20)
leg = legend('$\kappa = 0.25$','$\kappa = 0.275$','$\kappa = 0.30$','$\kappa = 0.325$','$\kappa = 0.35$',...
    '$\kappa = 0.9$','$\kappa = 0.95$','$\kappa = 1.0$','$\kappa = 1.05$','$\kappa = 1.1$')
set(leg,'fontsize',20,'interpreter','latex')
grid on
%axis([40 1000 -2.5 5])


%%

stringText = []; %zeros(10,200);

disp('=================================================================== ')
disp('MSE Values')
for int_counter = 1:size(mean_MSD_mat_LMS_comp,2)
    [MSE_theory, MSD_theory] = evaluate_NLMS_MSE_MSD(mu_vec(int_counter), D, U_f, diag(variance_vector_matrix(:,int_counter)));
    
    MSE_simu = mean( mean_error_mat_LMS_comp((end-numbCheckedUpdates):end,int_counter)); %10*log10( mean( mean_error_mat_LMS_comp((end-numbCheckedUpdates):end,int_counter)) );
    MSD_simu = mean( mean_MSD_mat_LMS_comp((end-numbCheckedUpdates):end,int_counter)) ;
    disp( '--------------------------------------------------------' )
    disp( [ 'mu = ' num2str(mu_vec(int_counter)) ] );
    disp( [ 'Simulation - MSE = ' num2str(MSE_simu) ' - MSD = ' num2str(MSD_simu) ] );
    disp( [ 'Theory - MSE = ' num2str(MSE_theory) ' - MSD = ' num2str(MSD_theory) ] );
    disp( [ 'Error - MSE = ' num2str((MSE_theory-MSE_simu)/MSE_theory) ' - MSD = ' num2str((MSD_theory-MSD_simu)/MSD_theory) ] );
    currentString = [ '$' num2str(mu_vec(int_counter)) '$ & $' num2str(MSE_simu) '$ & $' num2str(MSE_theory)  '$ & $' num2str(MSD_simu) '$ & $' num2str(MSD_theory) '$  \\ ' ]; % \hline' ];
    disp(currentString)
    
    stringText = [stringText ' ' currentString];
    
    % ----------------------------------------
    % ----------------------------------------
    if( int_counter <= length(mu_vec)/2 )
        idx = 2;
    else
        idx = 8;
    end
    % ----------------------------------------
    insertData = true;
    switch( mu_vec(int_counter) )
        case 0.05
            lineIndex = 1;
        case 0.10
            lineIndex = 2;
        case 0.25
            lineIndex = 3;
        case 0.50
            lineIndex = 4;
        otherwise
            insertData = false;
    end
    % ----------------------------------------
    %numbDecimalPlaces = 4;
    if(insertData == true)
        articleTable(lineIndex,idx) = round(MSE_theory*10^4)/(10^4);
        articleTable(lineIndex,idx+1) = round(MSE_simu*10^4)/(10^4); 
        articleTable(lineIndex,idx+2) = 100*((articleTable(lineIndex,idx)-articleTable(lineIndex,idx+1))/articleTable(lineIndex,idx)); 
        articleTable(lineIndex,idx+3) = round(MSD_theory*10^4)/(10^4);
        articleTable(lineIndex,idx+4) = round(MSD_simu*10^4)/(10^4);
        articleTable(lineIndex,idx+5) = 100*((articleTable(lineIndex,idx+3)-articleTable(lineIndex,idx+4))/articleTable(lineIndex,idx+3));
    end
    % ----------------------------------------
    % ----------------------------------------
    
end
disp('=================================================================== ')

%% ===================================================================== %%
% ======================================================================= %
% == II - GSP LMS Algorithm ============================================= %
% ======================================================================= %

load('02-simu-DS-LMS')

disp('=================================================================== ')
disp(['Ensemble of ' num2str(ensemble) ' runs'])
disp('Steady-state values taken considering the average of the last 2000 time instants.')
disp('=================================================================== ')
disp(['mu value = ' num2str(mu_vec(1))])
disp('=================================================================== ')

%% --------------------------------------------------------------------- %%
t = 1:length( mean_MSD_mat_LMS_comp(:,1) );
% ----------------------------------------------------------------------- %
% -- MSD (Mean Square Deviation) ---------------------------------------- %
figure
for int_counter = 1:length(mu_vec) 
    plot( t, 10*log10( mean_MSD_mat_LMS_comp(:,int_counter) ), 'LineWidth', 2 )
    hold on
end
ylabel('MSD [dB]','Interpreter','latex','fontsize',20)
xlabel('Number of Iterations $k$','Interpreter','latex','fontsize',20)
leg = legend('$\kappa = 0.25$','$\kappa = 0.275$','$\kappa = 0.30$','$\kappa = 0.325$','$\kappa = 0.35$',...
    '$\kappa = 0.9$','$\kappa = 0.95$','$\kappa = 1.0$','$\kappa = 1.05$','$\kappa = 1.1$')
set(leg,'fontsize',20,'interpreter','latex')
grid on
%axis([40 1000 -2.5 5])

%%
stringText = [];
disp('=================================================================== ')
for int_counter = 1:size(mean_MSD_mat_LMS_comp,2)
    [MSE_theory, MSD_theory, MSD_theory_Lor] = evaluate_LMS_MSE_MSD(mu_vec(int_counter), D, U_f, diag(variance_vector_matrix(:,int_counter)));
    
    MSE_simu = mean( mean_error_mat_LMS_comp((end-numbCheckedUpdates):end,int_counter)); %10*log10( mean( mean_error_mat_LMS_comp((end-numbCheckedUpdates):end,int_counter)) );
    MSD_simu = mean( mean_MSD_mat_LMS_comp((end-numbCheckedUpdates):end,int_counter)) ;
    disp( '--------------------------------------------------------' )
    disp( [ 'mu = ' num2str(mu_vec(int_counter)) ] );
    disp( [ 'Simulation - MSE = ' num2str(MSE_simu) ' - MSD = ' num2str(MSD_simu) ] );
    disp( [ 'Theory - MSE = ' num2str(MSE_theory) ' - MSD = ' num2str(MSD_theory) ' - MSD_Lorenzo = ' num2str(MSD_theory_Lor) ] );
    disp( [ 'Error - MSE = ' num2str((MSE_theory-MSE_simu)/MSE_theory) ' - MSD = ' num2str((MSD_theory-MSD_simu)/MSD_theory) ] );
    currentString = [ '$' num2str(mu_vec(int_counter)) '$ & $' num2str(MSE_simu) '$ & $' num2str(MSE_theory)  '$ & $' num2str(MSD_simu) '$ & $' num2str(MSD_theory) '$ & $' num2str(MSD_theory_Lor) '$  \\ ' ]; %\hline' ];
    disp(currentString)
    
    stringText = [stringText ' ' currentString];
    
    % ----------------------------------------
    % ----------------------------------------
    if( int_counter <= length(mu_vec)/2 )
        idx = 2;
    else
        idx = 8;
    end
    % ----------------------------------------
    insertData = true;
    switch( mu_vec(int_counter) )
        case 0.2
            lineIndex = 5;
        case 0.5
            lineIndex = 6;
        case 1.0
            lineIndex = 7;
        otherwise
            insertData = false;
    end
    % ----------------------------------------
    %numbDecimalPlaces = 4;
    if(insertData == true)
        articleTable(lineIndex,idx) = round(MSE_theory*10^4)/(10^4);
        articleTable(lineIndex,idx+1) = round(MSE_simu*10^4)/(10^4); 
        articleTable(lineIndex,idx+2) = 100*((articleTable(lineIndex,idx)-articleTable(lineIndex,idx+1))/articleTable(lineIndex,idx)); 
        articleTable(lineIndex,idx+3) = round(MSD_theory*10^4)/(10^4);
        articleTable(lineIndex,idx+4) = round(MSD_simu*10^4)/(10^4);
        articleTable(lineIndex,idx+5) = 100*((articleTable(lineIndex,idx+3)-articleTable(lineIndex,idx+4))/articleTable(lineIndex,idx+3));
    end
    % ----------------------------------------
    % ----------------------------------------
end
disp('=================================================================== ')

%% ===================================================================== %%
% ======================================================================= %
% == III - GSP RLS Algorithm ============================================ %
% ======================================================================= %

load('02-simu-DS-RLS')

disp('=================================================================== ')
disp(['Ensemble of ' num2str(ensemble) ' runs'])
disp('Steady-state values taken considering the average of the last 2000 time instants.')
disp('=================================================================== ')
%disp(['beta value = ' num2str(beta_vec(1))])
%disp('=================================================================== ')

%% --------------------------------------------------------------------- %%
% t = 1:length( mean_MSD_mat_RLS_comp(:,1) );
% % ----------------------------------------------------------------------- %
% % -- MSD (Mean Square Deviation) ---------------------------------------- %
% figure
% for int_counter = 1:length(beta_vec) 
%     plot( t, 10*log10( mean_MSD_mat_RLS_comp(:,int_counter) ), 'LineWidth', 2 )
%     hold on
% end
% ylabel('MSD [dB]','Interpreter','latex','fontsize',20)
% xlabel('Number of Iterations $k$','Interpreter','latex','fontsize',20)
% leg = legend('$\kappa = 0.25$','$\kappa = 0.275$','$\kappa = 0.30$','$\kappa = 0.325$','$\kappa = 0.35$',...
%     '$\kappa = 0.9$','$\kappa = 0.95$','$\kappa = 1.0$','$\kappa = 1.05$','$\kappa = 1.1$')
% set(leg,'fontsize',20,'interpreter','latex')
% grid on


%%
stringText = [];
disp('=================================================================== ')
for int_counter = 1:size(mean_MSD_mat_RLS_comp,2)
    [MSE_theory, MSD_theory] = evaluate_RLS_MSE_MSD(beta_vec(int_counter), D, U_f, diag(variance_vector_matrix(:,int_counter)));
    
    MSE_simu = mean( mean_error_mat_RLS_comp((end-numbCheckedUpdates):end,int_counter)); %10*log10( mean( mean_error_mat_RLS_comp((end-numbCheckedUpdates):end,int_counter)) );
    MSD_simu = mean( mean_MSD_mat_RLS_comp((end-numbCheckedUpdates):end,int_counter)) ;
    disp( '--------------------------------------------------------' )
    disp( [ 'beta = ' num2str(beta_vec(int_counter)) ] );
    disp( [ 'Simulation - MSE = ' num2str(MSE_simu) ' - MSD = ' num2str(MSD_simu) ] );
    disp( [ 'Theory - MSE = ' num2str(MSE_theory) ' - MSD = ' num2str(MSD_theory) ] );
    disp( [ 'Error - MSE = ' num2str((MSE_theory-MSE_simu)/MSE_theory) ' - MSD = ' num2str((MSD_theory-MSD_simu)/MSD_theory) ] );
    currentString = [ '$' num2str(beta_vec(int_counter)) '$ & $' num2str(MSE_simu) '$ & $' num2str(MSE_theory)  '$ & $' num2str(MSD_simu) '$ & $' num2str(MSD_theory) '$  \\ ' ]; %\hline' ];
    disp(currentString)
    
    stringText = [stringText ' ' currentString];
    
    % ----------------------------------------
    % ----------------------------------------
    if( int_counter <= length(beta_vec)/2 )
        idx = 2;
    else
        idx = 8;
    end
    % ----------------------------------------
    insertData = true;
    switch( beta_vec(int_counter) )
        case 0.95
            lineIndex = 8;
        case 0.9
            lineIndex = 9;
        case 0.75
            lineIndex = 10;
        otherwise
            insertData = false;
    end
    % ----------------------------------------
    %numbDecimalPlaces = 4;
    if(insertData == true)
        articleTable(lineIndex,idx) = round(MSE_theory*10^4)/(10^4);
        articleTable(lineIndex,idx+1) = round(MSE_simu*10^4)/(10^4); 
        articleTable(lineIndex,idx+2) = 100*((articleTable(lineIndex,idx)-articleTable(lineIndex,idx+1))/articleTable(lineIndex,idx)); 
        articleTable(lineIndex,idx+3) = round(MSD_theory*10^4)/(10^4);
        articleTable(lineIndex,idx+4) = round(MSD_simu*10^4)/(10^4);
        articleTable(lineIndex,idx+5) = 100*((articleTable(lineIndex,idx+3)-articleTable(lineIndex,idx+4))/articleTable(lineIndex,idx+3));
    end
    % ----------------------------------------
    % ----------------------------------------
end
disp('=================================================================== ')

%% ===================================================================== %%

table_string = [];
for lineCounter = 1:size(articleTable,1)
    line_string = [];
    for columnCounter = 1:14
        switch columnCounter
            case 1
                if (lineCounter == 1)
                    aux_string = [' \multirow{4}{*}{NLMS} ' ];
                    disp('Entrou')
                    aux_string
                else
                    if(lineCounter == 5)
                        aux_string = [' \multirow{3}{*}{LMS} ' ];
                        disp('Entrou')
                    aux_string
                    else
                        if(lineCounter == 8)
                            aux_string = [' \multirow{3}{*}{RLS} ' ];
                         else
                             aux_string = [' '];
                        end
                    end
                end
            case {2} %, 5, 7, 9, 11, 13, 15}
                aux_string = sprintf('%0.2f', articleTable(lineCounter,1) );
            case {5,8,11,14} 
                aux_string = sprintf('%0.3f', articleTable(lineCounter,columnCounter-1) );
            otherwise
                aux_string = sprintf('%0.4f', articleTable(lineCounter,columnCounter-1) );
        end
        aux_string
        if(columnCounter == 1)
            line_string = aux_string;
        else
            line_string = [line_string ' & ' aux_string];
        end
    end
    % ---------------------------------
    if ( (lineCounter == 4) || (lineCounter == 7) || (lineCounter == 10) )
        line_string = [line_string ' \\ \hline '];
    else
        line_string = [line_string ' \\ '];
    end
    % ---------------------------------   
    table_string = [ table_string line_string ];
end