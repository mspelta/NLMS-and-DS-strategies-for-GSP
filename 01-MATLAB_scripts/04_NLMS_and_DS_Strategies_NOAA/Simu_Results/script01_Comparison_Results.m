% ======================================================================= %
% == COPPE/UFRJ - Programa de Engenharia Eletrica (PEE) ================= %
% == Script: script01_Comparison_Rsults ================================= %
% == Responsible: Marcelo Jorge Mendes Spelta - Date: 2019/06/20 ======== %
% == E-mail: marcelo.spelta@smt.ufrj.br ================================= %
% ======================================================================= %

clc; clear;

% ======================================================================= %
% -- Loading original graph signal, bandlimited graph signal, U_f and  -- %
% -- D_s. Graph signal is obtained from the NOAA dataset. --------------- %
load('../General_Temperature_Data')
% ======================================================================= %

%% --------------------------------------------------------------------- %%
% -- Comparison between the LMS and NLMS algorithms
load('01-simu-StaticComparison_LMS-RLS-NLMS')

%% --------------------------------------------------------------------- %%
t = 1:95;
% ----------------------------------------------------------------------- %
% -- MSD (Mean Square Deviation) ---------------------------------------- %
figure
for int_counter = 1:3 %1:length(alg_param_vec) 
    sampleRate = 1;
    y_values = mean_MSD_mat_comp(1:sampleRate:length(t),int_counter);
    x_values = t(1:sampleRate:end);
    
    if(int_counter == 3)
        plot( x_values, 10*log10( y_values ), ':', 'LineWidth', 2 )
    else
        plot( x_values, 10*log10( y_values ), 'LineWidth', 2 )
    end
    hold on
end
ylabel('MSD$_{\rm G}$ [dB]','Interpreter','latex','fontsize',20)
xlabel('Number of Iterations $k$','Interpreter','latex','fontsize',20)
leg = legend('LMS','RLS','NLMS')
set(leg,'fontsize',18,'interpreter','latex')
axis([0 95 16.5 36]) 
grid on

%%
figProp = struct('size',22,'font','Times','lineWidth',4,'figDim',[1 1 800 350]);%600 400]);
figFileName = './Figures/StaticComparison_LMS_NLMS_RLS_MSD';
formatFig(gcf,figFileName,'en',figProp);

% ======================================================================= %
%% --------------------------------------------------------------------- %%
% -- Comparison between the LMS and NLMS algorithms
load('01-simu-DynamicComparison_LMS-RLS-NLMS')

%% --------------------------------------------------------------------- %%
% ----------------------------------------------------------------------- %
% -- MSD (Mean Square Deviation) ---------------------------------------- %
figure
x_values = 1:95;
sampled_position = 100; 
plot( x_values, completeDatasetMatrix(sampled_position,:) , 'k', 'LineWidth', 1 )
hold on
for int_counter = 1:3 
    sampleRate = 1;
    y_values = mean_x_vec_matrix(sampled_position,:,int_counter); 
    
    if(int_counter == 3)
        plot( x_values, y_values , ':', 'LineWidth', 2 )
    else
        plot( x_values, y_values , 'LineWidth', 2 )
    end
    hold on
end
ylabel('Temperature (°C)','Interpreter','latex','fontsize',20)
xlabel('Number of Iterations $k$','Interpreter','latex','fontsize',20)
leg = legend('Original','LMS','RLS','NLMS')
set(leg,'fontsize',18,'interpreter','latex')
axis([0 95 -2 +12]) 
grid on

%% --------------------------------------------------------------------- %%
% ----------------------------------------------------------------------- % 
figProp = struct('size',22,'font','Times','lineWidth',4,'figDim',[1 1 800 350]);
figFileName = './Figures/DynamicComparison_LMS_NLMS_RLS_MSD_100';
formatFig(gcf,figFileName,'en',figProp);

%% --------------------------------------------------------------------- %%
% ----------------------------------------------------------------------- %
% -- MSD (Mean Square Deviation) ---------------------------------------- %
figure
x_values = 1:95;
sampled_position = 119; 
plot( x_values, completeDatasetMatrix(sampled_position,:) , 'k', 'LineWidth', 1 )
hold on
for int_counter = 1:3 
    sampleRate = 1;
    y_values = mean_x_vec_matrix(sampled_position,:,int_counter); 
    
    if(int_counter == 3)
        plot( x_values, y_values , ':', 'LineWidth', 2 )
    else
        plot( x_values, y_values , 'LineWidth', 2 )
    end
    hold on
end
ylabel('Temperature (°C)','Interpreter','latex','fontsize',20)
xlabel('Number of Iterations $k$','Interpreter','latex','fontsize',20)
leg = legend('Original','LMS','RLS','NLMS')
set(leg,'fontsize',18,'interpreter','latex')
axis([0 95 +7 +24.5]) 
grid on

%% --------------------------------------------------------------------- %%
% ----------------------------------------------------------------------- % 
figProp = struct('size',22,'font','Times','lineWidth',4,'figDim',[1 1 800 350]);
figFileName = './Figures/DynamicComparison_LMS_NLMS_RLS_MSD_119';
formatFig(gcf,figFileName,'en',figProp);

% == END OF SCRIPT ====================================================== %
% ======================================================================= %