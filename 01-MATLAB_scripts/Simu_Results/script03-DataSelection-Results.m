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

%% ===================================================================== %%
% ----------------------------------------------
articleTable = zeros(6,15);
kappa_vec = [ 3 3.5 3.75 1 1.1 1.15];
for i = 1:6
    if( i <= 3)
        articleTable(i,1) = 1;
    else
        articleTable(i,1) = 2;
    end
    
    articleTable(i,2) = kappa_vec(i);
    
    articleTable(i,3) = evaluate_updateProbability(articleTable(i,1), kappa_vec(i), D)*100;
end
% ----------------------------------------------

%% ===================================================================== %%

load('03-simu-DS-NLMS_DataSel')

disp('=================================================================== ')
disp(['Ensemble of ' num2str(ensemble) ' runs'])
disp('Steady-state values taken considering the average of the last 2000 time instants.')
disp('=================================================================== ')
disp(['mu value = ' num2str(mu)])
disp('=================================================================== ')

%% --------------------------------------------------------------------- %%
t = 1:500; %length( mean_MSD_mat_LMS_comp(:,1) );
% ----------------------------------------------------------------------- %
% -- MSD (Mean Square Deviation) ---------------------------------------- %
figure
for int_counter = [ 1 3 5 6 ] %1:6 %1:(length(mu_vec)/2) 
    if( ( int_counter == 3 ) || ( int_counter == 6 ) )
        plot( t, 10*log10( mean_MSD_mat_LMS_comp(t,int_counter) ), 'LineWidth', 2 ) % '--', 'LineWidth', 2 )
    else
         plot( t, 10*log10( mean_MSD_mat_LMS_comp(t,int_counter) ), 'LineWidth', 2 )
    end
    hold on
end
ylabel('MSD [dB]','Interpreter','latex','fontsize',20)
xlabel('Number of Iterations $k$','Interpreter','latex','fontsize',20)
%leg = legend('$\kappa = 0.250$','$\kappa = 0.275$','$\kappa = 0.300$','$\kappa = 0.325$','$\kappa = 0.350$','$\kappa = 0.375$')
leg = legend('$\kappa = 2.50$','$\kappa = 3.00$','$\kappa = 3.50$','$\kappa = 3.75$')
set(leg,'fontsize',28,'interpreter','latex')
grid on

axis([40 500 -5 10])  %30])

%%
figProp = struct('size',24,'font','Times','lineWidth',4,'figDim',[1 1 800 400]);%600 400]);
figFileName = './Figures/DataSelection_NLMS_CW';
formatFig(gcf,figFileName,'en',figProp);

%% --------------------------------------------------------------------- %%
t = 1:400; %length( mean_MSD_mat_LMS_comp(:,1) );
% ----------------------------------------------------------------------- %
% -- MSD (Mean Square Deviation) ---------------------------------------- %
figure
for int_counter = [7 9 11 12] %7:12 %8:8 %12:12 %12 %(length(mu_vec)/2 + 1):(length(mu_vec)) 
    plot( t, 10*log10( mean_MSD_mat_LMS_comp(t,int_counter) ), 'LineWidth', 2 )
    hold on
end
ylabel('MSD [dB]','Interpreter','latex','fontsize',20)
xlabel('Number of Iterations $k$','Interpreter','latex','fontsize',20)
leg = legend('$\kappa = 0.90$','$\kappa = 1.00$','$\kappa = 1.10$','$\kappa = 1.15$')
set(leg,'fontsize',20,'interpreter','latex')
grid on

axis([55 400 -4.5 2])

%%
figProp = struct('size',24,'font','Times','lineWidth',4,'figDim',[1 1 800 400]);%600 400]);
figFileName = './Figures/DataSelection_NLMS_l2-Norm';
formatFig(gcf,figFileName,'en',figProp);

%%

stringText = []; %zeros(10,200);

disp('=================================================================== ')
%disp('MSE Values')
for int_counter = 1:size(mean_MSD_mat_LMS_comp,2)
    [ updateProbability ] = evaluate_updateProbability(DS_strategy_vec(int_counter), kappa_vec(int_counter), D);
    
    disp( '--------------------------------------------------------' )
    disp( [ 'DS Strategy ' num2str(DS_strategy_vec(int_counter)) ' - kappa = ' num2str(kappa_vec(int_counter)) ] );
    disp( [ 'Theory Update Rate[%] = ' num2str(updateProbability*100) ' -- Simulation Update Rate[%] = ' num2str(update_counter_vec_LMS_comp(int_counter)*100) ] );
    relError = (updateProbability - update_counter_vec_LMS_comp(int_counter))/updateProbability;
    disp( [ 'Error  = ' num2str(relError) ] );
    
    if(int_counter <= length(kappa_vec)/2 )
        cov_string = '(2)';
    else
        cov_string = '(4)';
    end
    
    if ( mod(int_counter,2) == 1)
        currentString = [ '(' num2str(DS_strategy_vec(int_counter)) ') & ' cov_string ' & ' num2str(kappa_vec(int_counter)) ...
            ' & ' num2str(update_counter_vec_LMS_comp(int_counter)*100) ' & ' num2str(updateProbability*100) ' & ']; % % \hline' ];
    else
        currentString = [  num2str(kappa_vec(int_counter)) ' & ' num2str(update_counter_vec_LMS_comp(int_counter)*100) ...
            ' & ' num2str(updateProbability*100) ' \\ ']; % % \hline' ];
    end
    %disp(currentString)
    
    stringText = [stringText ' ' currentString];
    
    % ----------------------------------------
    % ----------------------------------------
    if( cov_string == '(2)')
        idx = 4;
    else
        idx = 10;
    end
    % ----------------------------------------
    insertData = true;
    switch( kappa_vec(int_counter) )
        case 3
            lineIndex = 1;
        case 3.5
            lineIndex = 2;
        case 3.75
            lineIndex = 3;
        case 1.00
            lineIndex = 4;
        case 1.10
            lineIndex = 5;
        case 1.15
            lineIndex = 6;
        otherwise
            insertData = false;
    end
    % ----------------------------------------
    if(insertData == true)
        articleTable(lineIndex,idx) = round((update_counter_vec_LMS_comp(int_counter)*100)*10^3)/(10^3);
        shortUpdProb = round((updateProbability*100)*10^3)/(10^3);
        relError = (shortUpdProb - articleTable(lineIndex,idx))/shortUpdProb;
        articleTable(lineIndex,idx+1) = relError*100;
    end
    % ----------------------------------------
    % ----------------------------------------
    
    
    % ----------------------------------------
end
disp('=================================================================== ')

%% ===================================================================== %%
% ======================================================================= %
% == II - GSP LMS Algorithm ============================================= %
% ======================================================================= %

load('03-simu-DS-LMS')

disp('=================================================================== ')
disp(['Ensemble of ' num2str(ensemble) ' runs'])
disp('Steady-state values taken considering the average of the last 500 time instants.')
disp('=================================================================== ')
disp(['mu value = ' num2str(mu)])
disp('=================================================================== ')
%% --------------------------------------------------------------------- %%
t = 1:length( mean_MSD_mat_LMS_comp(:,1) );
% ----------------------------------------------------------------------- %
% -- MSD (Mean Square Deviation) ---------------------------------------- %
figure
for int_counter = 15:16 %1:length(mu_vec) 
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
    
    [ updateProbability ] = evaluate_updateProbability(DS_strategy_vec(int_counter), kappa_vec(int_counter), D);
    
    disp( '--------------------------------------------------------' )
    disp( [ 'DS Strategy ' num2str(DS_strategy_vec(int_counter)) ' - kappa = ' num2str(kappa_vec(int_counter)) ] );
    disp( [ 'Theory Update Rate[%] = ' num2str(updateProbability*100) ' -- Simulation Update Rate[%] = ' num2str(update_counter_vec_LMS_comp(int_counter)*100) ] );
    relError = (updateProbability - update_counter_vec_LMS_comp(int_counter))/updateProbability;
    disp( [ 'Error  = ' num2str(relError) ] );
    
    if(int_counter <= length(kappa_vec)/2 )
        cov_string = '(2)';
    else
        cov_string = '(4)';
    end
    
    if ( mod(int_counter,2) == 1)
        currentString = [ '(' num2str(DS_strategy_vec(int_counter)) ') & ' cov_string ' & ' num2str(kappa_vec(int_counter)) ...
            ' & ' num2str(update_counter_vec_LMS_comp(int_counter)*100) ' & ' num2str(updateProbability*100) ' & ']; 
    else
        currentString = [  num2str(kappa_vec(int_counter)) ' & ' num2str(update_counter_vec_LMS_comp(int_counter)*100) ...
            ' & ' num2str(updateProbability*100) ' \\ ']; % % \hline' ];
    end
   
    %disp(currentString)
    
    stringText = [stringText ' ' currentString];

    % ----------------------------------------
    % ----------------------------------------
    if( cov_string == '(2)')
        idx = 6;
    else
        idx = 12;
    end
    % ----------------------------------------
    insertData = true;
    switch( kappa_vec(int_counter) )
        case 3
            lineIndex = 1;
        case 3.5
            lineIndex = 2;
        case 3.75
            lineIndex = 3;
        case 1.00
            lineIndex = 4;
        case 1.10
            lineIndex = 5;
        case 1.15
            lineIndex = 6;
        otherwise
            insertData = false;
    end
    % ----------------------------------------
    if(insertData == true)
        articleTable(lineIndex,idx) = round((update_counter_vec_LMS_comp(int_counter)*100)*10^3)/(10^3);
        shortUpdProb = round((updateProbability*100)*10^3)/(10^3);
        relError = (shortUpdProb - articleTable(lineIndex,idx))/shortUpdProb;
        articleTable(lineIndex,idx+1) = relError*100;
    end
    % ----------------------------------------
    % ----------------------------------------
    
end
disp('=================================================================== ')

%% ===================================================================== %%
% ======================================================================= %
% == III - GSP RLS Algorithm ============================================ %
% ======================================================================= %

load('03-simu-DS-RLS')

disp('=================================================================== ')
disp(['Ensemble of ' num2str(ensemble) ' runs'])
disp('Steady-state values taken considering the average of the last 500 time instants.')
disp('=================================================================== ')
%disp(['beta value = ' num2str(beta_vec(1))])
%disp('=================================================================== ')
%% --------------------------------------------------------------------- %%
t = 1:length( mean_MSD_mat_RLS_comp(:,1) );
% ----------------------------------------------------------------------- %
% -- MSD (Mean Square Deviation) ---------------------------------------- %
figure
for int_counter = 1:length(beta_vec) 
    plot( t, 10*log10( mean_MSD_mat_RLS_comp(:,int_counter) ), 'LineWidth', 2 )
    hold on
end
ylabel('MSD [dB]','Interpreter','latex','fontsize',20)
xlabel('Number of Iterations $k$','Interpreter','latex','fontsize',20)
leg = legend('$\kappa = 0.25$','$\kappa = 0.275$','$\kappa = 0.30$','$\kappa = 0.325$','$\kappa = 0.35$',...
    '$\kappa = 0.9$','$\kappa = 0.95$','$\kappa = 1.0$','$\kappa = 1.05$','$\kappa = 1.1$')
set(leg,'fontsize',20,'interpreter','latex')
grid on


%%
stringText = [];
disp('=================================================================== ')
for int_counter = 1:size(mean_MSD_mat_RLS_comp,2)
    
    [ updateProbability ] = evaluate_updateProbability(DS_strategy_vec(int_counter), kappa_vec(int_counter), D);
    
    disp( '--------------------------------------------------------' )
    disp( [ 'DS Strategy ' num2str(DS_strategy_vec(int_counter)) ' - kappa = ' num2str(kappa_vec(int_counter)) ] );
    disp( [ 'Theory Update Rate[%] = ' num2str(updateProbability*100) ' -- Simulation Update Rate[%] = ' num2str(update_counter_vec_RLS_comp(int_counter)*100) ] );
    relError = (updateProbability - update_counter_vec_RLS_comp(int_counter))/updateProbability;
    disp( [ 'Error  = ' num2str(relError) ] );
    
    if(int_counter <= length(kappa_vec)/2 )
        cov_string = '(2)';
    else
        cov_string = '(4)';
    end
    
    if ( mod(int_counter,2) == 1)
        currentString = [ '(' num2str(DS_strategy_vec(int_counter)) ') & ' cov_string ' & ' num2str(kappa_vec(int_counter)) ...
            ' & ' num2str(update_counter_vec_RLS_comp(int_counter)*100) ' & ' num2str(updateProbability*100) ' & ']; 
    else
        currentString = [  num2str(kappa_vec(int_counter)) ' & ' num2str(update_counter_vec_RLS_comp(int_counter)*100) ...
            ' & ' num2str(updateProbability*100) ' \\ ']; % % \hline' ];
    end
    
    stringText = [stringText ' ' currentString];
    
    % ----------------------------------------
    % ----------------------------------------
    if( cov_string == '(2)')
        idx = 8;
    else
        idx = 14;
    end
    % ----------------------------------------
    insertData = true;
    switch( kappa_vec(int_counter) )
        case 3
            lineIndex = 1;
        case 3.5
            lineIndex = 2;
        case 3.75
            lineIndex = 3;
        case 1.00
            lineIndex = 4;
        case 1.10
            lineIndex = 5;
        case 1.15
            lineIndex = 6;
        otherwise
            insertData = false;
    end
    % ----------------------------------------
    if(insertData == true)
        articleTable(lineIndex,idx) = round((update_counter_vec_RLS_comp(int_counter)*100)*10^3)/(10^3);
        shortUpdProb = round((updateProbability*100)*10^3)/(10^3);
        relError = (shortUpdProb - articleTable(lineIndex,idx))/shortUpdProb;
        articleTable(lineIndex,idx+1) = relError*100;
    end
    % ----------------------------------------
    % ----------------------------------------
end
disp('=================================================================== ')

%% ===================================================================== %%

table_string = [];
for lineCounter = 1:size(articleTable,1)
    line_string = [];
    for columnCounter = 1:15
        switch columnCounter
            case 1
                aux_string = ['(' num2str(articleTable(lineCounter,columnCounter)) ')' ];
            case {2} %, 5, 7, 9, 11, 13, 15}
                aux_string = sprintf('%0.2f', articleTable(lineCounter,columnCounter));
            %case 5, 7, 9, 11, 13, 15
            otherwise
                aux_string = sprintf('%0.3f', articleTable(lineCounter,columnCounter));
        end
        if(columnCounter == 1)
            line_string = aux_string;
        else
            line_string = [line_string ' & ' aux_string];
        end
    end
    line_string = [line_string ' \\ '];
    
    table_string = [ table_string line_string ];
end