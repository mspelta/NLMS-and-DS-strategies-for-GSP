% ======================================================================= %
% == COPPE/UFRJ - Programa de Engenharia Eletrica (PEE) ================= %
% == Script: script02_Analysis_Results ================================== %
% == Responsible: Marcelo Jorge Mendes Spelta - Date: 2019/03/29 ======== %
% == E-mail: marcelo.spelta@smt.ufrj.br ================================= %
% ======================================================================= %

clc; clear;

numbCheckedUpdates = 1000;

% ======================================================================= %
% -- Loading original graph signal, bandlimited graph signal, U_f and  -- %
% -- D_s. Graph signal is obtained from the INMET dataset. -------------- %
load('../General_Bandlimited_GS_Data')
% ======================================================================= %

%% ===================================================================== %%
% == I - GSP NLMS Algorithm ============================================= %
% ======================================================================= %

load('02-simu-DS-NLMS')

% ----------------------------------------------
tableScenario2 = zeros(10,7);
tableScenario3 = zeros(10,7);
factors_vec = mu_vec(1:4);
tableScenario2(1:4,1) = factors_vec;
tableScenario3(1:4,1) = factors_vec;
% ----------------------------------------------

disp('=================================================================== ')
disp(['Ensemble of ' num2str(ensemble) ' runs'])
disp('Steady-state values taken considering the average of the last 2000 time instants.')
disp('=================================================================== ')

%% ===================================================================== %%
% == Visualization of the MSD curve for the NLMS algorithm. ============= %
% == (Illustration only) ================================================ %
% ----------------------------------------------------------------------- %

t = 1:length( mean_MSD_mat_LMS_comp(:,1) );
figure
for int_counter = 1:length(mu_vec) 
    plot( t, 10*log10( mean_MSD_mat_LMS_comp(:,int_counter) ), 'LineWidth', 2 )
    hold on
end
ylabel('MSD [dB]','Interpreter','latex','fontsize',20)
xlabel('Number of Iterations $k$','Interpreter','latex','fontsize',20)
leg = legend(['$\mu_{\rm N} =$' num2str(mu_vec(1)) ' - $(ii)$'],['$\mu_{\rm N} =$' num2str(mu_vec(2)) ' - $(ii)$'],['$\mu_{\rm N} =$' num2str(mu_vec(3)) ' - $(ii)$'],['$\mu_{\rm N} =$' num2str(mu_vec(4)) ' - $(ii)$'],...
    ['$\mu_{\rm N} =$' num2str(mu_vec(5)) ' - $(iii)$'],['$\mu_{\rm N} =$' num2str(mu_vec(6)) ' - $(iii)$'],['$\mu_{\rm N} =$' num2str(mu_vec(7)) ' - $(iii)$'],['$\mu_{\rm N} =$' num2str(mu_vec(8)) ' - $(iii)$'])
set(leg,'fontsize',20,'interpreter','latex')
grid on

%% ===================================================================== %%
% == Obtaining estimates of the stationary MSE and MSD for the GSP ====== %
% == NLMS algorithm. ==================================================== %

x_diff = graph_signal - bandlimited_graph;

disp('=================================================================== ')
disp('MSE Values')
for int_counter = 1:size(mean_MSD_mat_LMS_comp,2)
    [MSE_theory, MSD_theory] = evaluate_NLMS_MSE_MSD(mu_vec(int_counter), D_s, U_f, diag(variance_vector_matrix(:,int_counter)), x_diff);
    
    MSE_simu = mean( mean_error_mat_LMS_comp((end-numbCheckedUpdates):end,int_counter)); %10*log10( mean( mean_error_mat_LMS_comp((end-numbCheckedUpdates):end,int_counter)) );
    MSD_simu = mean( mean_MSD_mat_LMS_comp((end-numbCheckedUpdates):end,int_counter)) ;
    disp( '--------------------------------------------------------' )
    disp( [ 'mu = ' num2str(mu_vec(int_counter)) ] );
    disp( [ 'Simulation - MSE = ' num2str(MSE_simu) ' - MSD = ' num2str(MSD_simu) ] );
    disp( [ 'Theory - MSE = ' num2str(MSE_theory) ' - MSD = ' num2str(MSD_theory) ] );
    disp( [ 'Error - MSE = ' num2str(100*(MSE_theory-MSE_simu)/MSE_simu) '% - MSD = ' num2str(100*(MSD_theory-MSD_simu)/MSD_simu) '%' ] );
    
    % ------------------------------------------------------------------- %
    % -- Writing data to the table that will be used for summarizing results 
    if (int_counter <= length(mu_vec)/2)
        lineIndex = int_counter;
        tableScenario2(lineIndex,2) = round(MSE_theory*10^4)/(10^4);
        tableScenario2(lineIndex,3) = round(MSE_simu*10^4)/(10^4); 
        tableScenario2(lineIndex,4) = 100*((tableScenario2(lineIndex,2)-tableScenario2(lineIndex,3))/tableScenario2(lineIndex,3)); 
        tableScenario2(lineIndex,5) = round(MSD_theory*10^4)/(10^4);
        tableScenario2(lineIndex,6) = round(MSD_simu*10^4)/(10^4);
        tableScenario2(lineIndex,7) = 100*((tableScenario2(lineIndex,5)-tableScenario2(lineIndex,6))/tableScenario2(lineIndex,6));
    else
        lineIndex = int_counter - 4;
        tableScenario3(lineIndex,2) = round(MSE_theory*10^4)/(10^4);
        tableScenario3(lineIndex,3) = round(MSE_simu*10^4)/(10^4); 
        tableScenario3(lineIndex,4) = 100*((tableScenario3(lineIndex,2)-tableScenario3(lineIndex,3))/tableScenario3(lineIndex,3)); 
        tableScenario3(lineIndex,5) = round(MSD_theory*10^4)/(10^4);
        tableScenario3(lineIndex,6) = round(MSD_simu*10^4)/(10^4);
        tableScenario3(lineIndex,7) = 100*((tableScenario3(lineIndex,5)-tableScenario3(lineIndex,6))/tableScenario3(lineIndex,6));
    end
    
end
disp('=================================================================== ')

%% ===================================================================== %%
% == II - GSP LMS Algorithm ============================================= %
% ======================================================================= %

load('02-simu-DS-LMS')

tableScenario2(5:7,1) = mu_vec(1:3);
tableScenario3(5:7,1) = mu_vec(1:3);

disp('=================================================================== ')
disp(['Ensemble of ' num2str(ensemble) ' runs'])
disp('Steady-state values taken considering the average of the last 2000 time instants.')
disp('=================================================================== ')

%% ===================================================================== %%
% == Visualization of the MSD curve for the LMS algorithm. ============== %
% == (Illustration only) ================================================ %
% ----------------------------------------------------------------------- %

t = 1:length( mean_MSD_mat_LMS_comp(:,1) );

figure
for int_counter = 1:length(mu_vec) 
    plot( t, 10*log10( mean_MSD_mat_LMS_comp(:,int_counter) ), 'LineWidth', 2 )
    hold on
end
ylabel('MSD [dB]','Interpreter','latex','fontsize',20)
xlabel('Number of Iterations $k$','Interpreter','latex','fontsize',20)
leg = legend(['$\mu_{\rm L} =$' num2str(mu_vec(1)) ' - $(ii)$'],['$\mu_{\rm L} =$' num2str(mu_vec(2)) ' - $(ii)$'],['$\mu_{\rm L} =$' num2str(mu_vec(3)) ' - $(ii)$'],...
    ['$\mu_{\rm L} =$' num2str(mu_vec(4)) ' - $(iii)$'],['$\mu_{\rm L} =$' num2str(mu_vec(5)) ' - $(iii)$'],['$\mu_{\rm L} =$' num2str(mu_vec(6)) ' - $(iii)$'])
set(leg,'fontsize',20,'interpreter','latex')
grid on

%% ===================================================================== %%
% == Obtaining estimates of the stationary MSE and MSD for the GSP ====== %
% == LMS algorithm. ===================================================== %

disp('=================================================================== ')
for int_counter = 1:size(mean_MSD_mat_LMS_comp,2)
    [MSE_theory, MSD_theory, MSD_theory_Lor] = evaluate_LMS_MSE_MSD(mu_vec(int_counter), D_s, U_f, diag(variance_vector_matrix(:,int_counter)));
    
    MSE_simu = mean( mean_error_mat_LMS_comp((end-numbCheckedUpdates):end,int_counter)); %10*log10( mean( mean_error_mat_LMS_comp((end-numbCheckedUpdates):end,int_counter)) );
    MSD_simu = mean( mean_MSD_mat_LMS_comp((end-numbCheckedUpdates):end,int_counter)) ;
    disp( '--------------------------------------------------------' )
    disp( [ 'mu = ' num2str(mu_vec(int_counter)) ] );
    disp( [ 'Simulation - MSE = ' num2str(MSE_simu) ' - MSD = ' num2str(MSD_simu) ] );
    disp( [ 'Theory - MSE = ' num2str(MSE_theory) ' - MSD = ' num2str(MSD_theory) ' - MSD_Lorenzo = ' num2str(MSD_theory_Lor) ] );
    disp( [ 'Error - MSE = ' num2str(100*(MSE_theory-MSE_simu)/MSE_simu) '% - MSD = ' num2str(100*(MSD_theory-MSD_simu)/MSD_simu) '%' ] );
      
    % ------------------------------------------------------------------- %
    % -- Writing data to the table that will be used for summarizing results 
    if (int_counter <= length(mu_vec)/2)
        lineIndex = int_counter + 4;
        tableScenario2(lineIndex,2) = round(MSE_theory*10^4)/(10^4);
        tableScenario2(lineIndex,3) = round(MSE_simu*10^4)/(10^4); 
        tableScenario2(lineIndex,4) = 100*((tableScenario2(lineIndex,2)-tableScenario2(lineIndex,3))/tableScenario2(lineIndex,3)); 
        tableScenario2(lineIndex,5) = round(MSD_theory*10^4)/(10^4);
        tableScenario2(lineIndex,6) = round(MSD_simu*10^4)/(10^4);
        tableScenario2(lineIndex,7) = 100*((tableScenario2(lineIndex,5)-tableScenario2(lineIndex,6))/tableScenario2(lineIndex,6));
    else
        lineIndex = int_counter + 1;
        tableScenario3(lineIndex,2) = round(MSE_theory*10^4)/(10^4);
        tableScenario3(lineIndex,3) = round(MSE_simu*10^4)/(10^4); 
        tableScenario3(lineIndex,4) = 100*((tableScenario3(lineIndex,2)-tableScenario3(lineIndex,3))/tableScenario3(lineIndex,3)); 
        tableScenario3(lineIndex,5) = round(MSD_theory*10^4)/(10^4);
        tableScenario3(lineIndex,6) = round(MSD_simu*10^4)/(10^4);
        tableScenario3(lineIndex,7) = 100*((tableScenario3(lineIndex,5)-tableScenario3(lineIndex,6))/tableScenario3(lineIndex,6));
    end
    
    % ----------------------------------------
    % ----------------------------------------
end
disp('=================================================================== ')

%% ===================================================================== %%
% == III - GSP RLS Algorithm ============================================ %
% ======================================================================= %

load('02-simu-DS-RLS')

tableScenario2(8:10,1) = beta_vec(1:3);
tableScenario3(8:10,1) = beta_vec(1:3);

disp('=================================================================== ')
disp(['Ensemble of ' num2str(ensemble) ' runs'])
disp('Steady-state values taken considering the average of the last 2000 time instants.')
disp('=================================================================== ')


%% ===================================================================== %%
% == Visualization of the MSD curve for the RLS algorithm. ============== %
% == (Illustration only) ================================================ %
% ----------------------------------------------------------------------- %

t = 1:length( mean_MSD_mat_RLS_comp(:,1) );

figure
for int_counter = 1:length(beta_vec) 
    plot( t, 10*log10( mean_MSD_mat_RLS_comp(:,int_counter) ), 'LineWidth', 2 )
    hold on
end
ylabel('MSD [dB]','Interpreter','latex','fontsize',20)
xlabel('Number of Iterations $k$','Interpreter','latex','fontsize',20)
leg = legend(['$\beta_{\rm R} =$' num2str(beta_vec(1)) ' - $(ii)$'],['$\beta_{\rm R} =$' num2str(beta_vec(2)) ' - $(ii)$'],['$\beta_{\rm R} =$' num2str(beta_vec(3)) ' - $(ii)$'],...
    ['$\beta_{\rm R} =$' num2str(beta_vec(4)) ' - $(iii)$'],['$\beta_{\rm R} =$' num2str(beta_vec(5)) ' - $(iii)$'],['$\beta_{\rm R} =$' num2str(beta_vec(6)) ' - $(iii)$'])
set(leg,'fontsize',20,'interpreter','latex')
grid on


%% ===================================================================== %%
% == Obtaining estimates of the stationary MSE and MSD for the GSP ====== %
% == NLMS algorithm. ==================================================== %

stringText = [];
disp('=================================================================== ')
for int_counter = 1:size(mean_MSD_mat_RLS_comp,2)
    [MSE_theory, MSD_theory] = evaluate_RLS_MSE_MSD(beta_vec(int_counter), D_s, U_f, diag(variance_vector_matrix(:,int_counter)));
    
    MSE_simu = mean( mean_error_mat_RLS_comp((end-numbCheckedUpdates):end,int_counter)); %10*log10( mean( mean_error_mat_RLS_comp((end-numbCheckedUpdates):end,int_counter)) );
    MSD_simu = mean( mean_MSD_mat_RLS_comp((end-numbCheckedUpdates):end,int_counter)) ;
    disp( '--------------------------------------------------------' )
    disp( [ 'beta = ' num2str(beta_vec(int_counter)) ] );
    disp( [ 'Simulation - MSE = ' num2str(MSE_simu) ' - MSD = ' num2str(MSD_simu) ] );
    disp( [ 'Theory - MSE = ' num2str(MSE_theory) ' - MSD = ' num2str(MSD_theory) ] );
    disp( [ 'Error - MSE = ' num2str(100*(MSE_theory-MSE_simu)/MSE_simu) '% - MSD = ' num2str(100*(MSD_theory-MSD_simu)/MSD_simu) '%'] );
    
    % ------------------------------------------------------------------- %
    % -- Writing data to the table that will be used for summarizing results 
    if (int_counter <= length(beta_vec)/2)
        lineIndex = int_counter + 7;
        tableScenario2(lineIndex,2) = round(MSE_theory*10^4)/(10^4);
        tableScenario2(lineIndex,3) = round(MSE_simu*10^4)/(10^4); 
        tableScenario2(lineIndex,4) = 100*((tableScenario2(lineIndex,2)-tableScenario2(lineIndex,3))/tableScenario2(lineIndex,3)); 
        tableScenario2(lineIndex,5) = round(MSD_theory*10^4)/(10^4);
        tableScenario2(lineIndex,6) = round(MSD_simu*10^4)/(10^4);
        tableScenario2(lineIndex,7) = 100*((tableScenario2(lineIndex,5)-tableScenario2(lineIndex,6))/tableScenario2(lineIndex,6));
    else
        lineIndex = int_counter + 4;
        tableScenario3(lineIndex,2) = round(MSE_theory*10^4)/(10^4);
        tableScenario3(lineIndex,3) = round(MSE_simu*10^4)/(10^4); 
        tableScenario3(lineIndex,4) = 100*((tableScenario3(lineIndex,2)-tableScenario3(lineIndex,3))/tableScenario3(lineIndex,3)); 
        tableScenario3(lineIndex,5) = round(MSD_theory*10^4)/(10^4);
        tableScenario3(lineIndex,6) = round(MSD_simu*10^4)/(10^4);
        tableScenario3(lineIndex,7) = 100*((tableScenario3(lineIndex,5)-tableScenario3(lineIndex,6))/tableScenario3(lineIndex,6));
    end
    
end
disp('=================================================================== ')

%% ===================================================================== %%
% == Editing data results in order to display them in a proper form for = %
% == the MSc dissertation. ============================================== %

table2 = []; table3 = [];

for ext_counter = 2:3 % 2 for tableScenario2 and 3 for tableScenario3
    for lineCounter = 1:size(tableScenario2,1)
        line_string = [];
        
        if (ext_counter == 2)
            lineValues = tableScenario2(lineCounter,:);
        else
            lineValues = tableScenario3(lineCounter,:);
        end
        
        for columnCounter = 1:8
            switch columnCounter
                case 1
                    if (lineCounter == 1)
                        aux_string = [' \multirow{4}{*}{NLMS} ' ];
                    else
                        if(lineCounter == 5)
                            aux_string = [' \multirow{3}{*}{LMS} ' ];
                        else
                            if(lineCounter == 8)
                                aux_string = [' \multirow{3}{*}{RLS} ' ];
                             else
                                 aux_string = [' '];
                            end
                        end
                    end
                case {2} % Algorithm factors
                    aux_string = sprintf('%0.2f', lineValues(1) );
                case {5,8} 
                    aux_string = sprintf('%0.3f', lineValues(columnCounter-1) );
                otherwise
                    aux_string = sprintf('%0.4f', lineValues(columnCounter-1) );
            end
            if(columnCounter == 1)
                line_string = aux_string;
            else
                line_string = [line_string ' & ' aux_string];
            end
        end
        % ---------------------------------
        if ( (lineCounter == 4) || (lineCounter == 7)  )
            line_string = [line_string ' \\ \hline '];
        else
            line_string = [line_string ' \\ '];
        end
        % ---------------------------------
        if (ext_counter == 2)
            table2 = [ table2 line_string ];
        else
            table3 = [ table3 line_string ];
        end
    end
end

% == END OF SCRIPT ====================================================== %
% ======================================================================= %