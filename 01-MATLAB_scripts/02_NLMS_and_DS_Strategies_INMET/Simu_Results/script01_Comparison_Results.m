%% ===================================================================== %%
% ======================================================================= %
% == I - GSP NLMS Algorithm ============================================= %
% ======================================================================= %

clc; clear;

numbCheckedUpdates = 1000;

% ======================================================================= %
% -- Loading original graph signal, bandlimited graph signal, U_f and  -- %
% -- D_s. Graph signal is obtained from the INMET dataset. -------------- %
load('../General_Bandlimited_GS_Data')
% ======================================================================= %


%% --------------------------------------------------------------------- %%
% -- Comparison between the LMS and NLMS algorithms
load('01-simu-Comparison_LMS-RLS-NLMS')

%% --------------------------------------------------------------------- %%
t = 1:4500; %length( mean_MSD_mat_comp(:,1) );
% ----------------------------------------------------------------------- %
% -- MSD (Mean Square Deviation) ---------------------------------------- %
figure
for int_counter = [ 1 3 4 6 7 9] %1:6 %length(alg_param_vec) 
%     switch(int_counter)
%         case 1
%             lineDetail = '-';
%         case 2 
%             lineDetail = '--';
%         case 3 
%             lineDetail = ':';
%         case 4
%             lineDetail = '-.';
%         case 5 
%             lineDetail = 's';
%         case 6 
%             lineDetail = 'd';
%     end
    sampleRate = 1;
    y_values = mean_MSD_mat_comp(1:sampleRate:length(t),int_counter);
    x_values = t(1:sampleRate:end);
    
    if(( int_counter == 4)||( int_counter == 6))
        plot( x_values, 10*log10( y_values ), ':', 'LineWidth', 2 )
    else
        if(( int_counter == 7)||( int_counter == 9))
            plot( x_values, 10*log10( y_values ), ':', 'LineWidth', 2 )
        else
            plot( x_values, 10*log10( y_values ), 'LineWidth', 2 )
        end
    end
    %plot( x_values, 10*log10( y_values ), lineDetail, 'LineWidth', 2 )
    hold on
end
ylabel('MSD$_{\rm G}$ [dB]','Interpreter','latex','fontsize',20)
xlabel('Number of Iterations $k$','Interpreter','latex','fontsize',20)
leg = legend('(L.I)','(N.I)','(L.II)','(N.II)','(L.III)','(N.III)')
set(leg,'fontsize',18,'interpreter','latex','NumColumns',2)
grid on
axis([0 4000 -20 52.5])

%%
figProp = struct('size',22,'font','Times','lineWidth',4,'figDim',[1 1 800 350]);%600 400]);
figFileName = './Figures/Comparison_LMS_NLMS_MSD';
formatFig(gcf,figFileName,'en',figProp);

%% --------------------------------------------------------------------- %%
t = 1:4500; %length( mean_MSD_mat_comp(:,1) );
% ----------------------------------------------------------------------- %
% -- MSD (Mean Square Deviation) ---------------------------------------- %
figure
aux_int_counter = 1;
for int_counter = [ 1 3 4 6 7 9] %1:6 %length(alg_param_vec)
    subplot(2,2,floor((aux_int_counter+1)/2))
%     switch(int_counter)
%         case 1
%             lineDetail = '-';
%         case 2 
%             lineDetail = '--';
%         case 3 
%             lineDetail = ':';
%         case 4
%             lineDetail = '-.';
%         case 5 
%             lineDetail = 's';
%         case 6 
%             lineDetail = 'd';
%     end
    sampleRate = 1;
    %y_values = mean_error_mat_comp(1:sampleRate:length(t),int_counter); %mean_orig_MSD_mat_comp(1:sampleRate:length(t),int_counter);
    y_values = mean_orig_MSD_mat_comp(1:sampleRate:length(t),int_counter);
    x_values = t(1:sampleRate:end);
    
    if(( int_counter == 4)||( int_counter == 6))
        plot( x_values, 10*log10( y_values ), ':', 'LineWidth', 2 )
    else
        if(( int_counter == 7)||( int_counter == 9))
            plot( x_values, 10*log10( y_values ), ':', 'LineWidth', 2 )
        else
            plot( x_values, 10*log10( y_values ), 'LineWidth', 2 )
        end
    end
    %plot( x_values, 10*log10( y_values ), lineDetail, 'LineWidth', 2 )
    hold on
    
    aux_int_counter = aux_int_counter + 1;
end
ylabel('MSD$_{\rm G}$ [dB]','Interpreter','latex','fontsize',20)
xlabel('Number of Iterations $k$','Interpreter','latex','fontsize',20)
leg = legend('(L.I)','(N.I)','(L.II)','(N.II)','(L.III)','(N.III)')
set(leg,'fontsize',18,'interpreter','latex','NumColumns',2)
grid on
axis([0 4000 10 52.5])

%%
figProp = struct('size',22,'font','Times','lineWidth',4,'figDim',[1 1 800 350]);%600 400]);
figFileName = './Figures/Comparison_LMS_NLMS_orig_MSD';
formatFig(gcf,figFileName,'en',figProp);

%% --------------------------------------------------------------------- %%
t = 1:4500; %length( mean_MSD_mat_comp(:,1) );
% ----------------------------------------------------------------------- %
% -- MSD (Mean Square Deviation) ---------------------------------------- %
figure
for int_counter = [4 5 6] % 1:length(alg_param_vec) 
%     switch(int_counter)
%         case 1
%             lineDetail = '-';
%         case 2 
%             lineDetail = '--';
%         case 3 
%             lineDetail = ':';
%         case 4
%             lineDetail = '-.';
%         case 5 
%             lineDetail = 's';
%         case 6 
%             lineDetail = 'd';
%     end
    sampleRate = 1;
    y_values = mean_MSD_mat_comp(1:sampleRate:length(t),int_counter);
    x_values = t(1:sampleRate:end);
    
    if(int_counter == 6)
        plot( x_values, 10*log10( y_values ), ':', 'LineWidth', 2 )
    else
        plot( x_values, 10*log10( y_values ), 'LineWidth', 2 )
    end
    %plot( x_values, 10*log10( y_values ), lineDetail, 'LineWidth', 2 )
    hold on
end
ylabel('MSD$_{\rm G}$ [dB]','Interpreter','latex','fontsize',20)
xlabel('Number of Iterations $k$','Interpreter','latex','fontsize',20)
leg = legend('(L.II)','(R.II)','(N.II)')
set(leg,'fontsize',18,'interpreter','latex')
axis([0 3500 24 45]) %([0 4000 -7.5 52.5])
grid on

%%
figProp = struct('size',22,'font','Times','lineWidth',4,'figDim',[1 1 800 350]);%600 400]);
figFileName = './Figures/Comparison_LMS_NLMS_RLS_MSD';
formatFig(gcf,figFileName,'en',figProp);

%%
figure
plot( mean_elapsedTime_mat_comp(:,1:3) ) %3:3:9) ) 

%% ===================================================================== %%
% ======================================================================= %
% == Evaluating the Average Stationary MSD and the Average Amount of ==== %
% == Iterations until the Algorithm Convergence. ======================== %
% ======================================================================= %

t = 1:length( mean_MSD_mat_comp(:,1) );

conv_time_vec = zeros(1,9);
conv_flags = zeros(1,9);

SteadyState_MSD_vec = zeros(1,9);

figure

for ext_loop = 0:2

    SteadyState_MSD_vec(ext_loop*3 + 1) =  mean( mean_MSD_mat_comp((end-numbCheckedUpdates):end, ext_loop*3 + 1) );  
    SteadyState_MSD_vec(ext_loop*3 + 2) =  mean( mean_MSD_mat_comp((end-numbCheckedUpdates):end, ext_loop*3 + 2) ); 
    SteadyState_MSD_vec(ext_loop*3 + 3) =  mean( mean_MSD_mat_comp((end-numbCheckedUpdates):end, ext_loop*3 + 3) ); 
    
    SteadyState_MSD = mean( SteadyState_MSD_vec( (ext_loop*3 + 1):(ext_loop*3 + 3) ) )
    SS_MSD_threshold = SteadyState_MSD*1.025; % Threshold value is 2.5% above the stationaty MSD

    for i = 1:length(t)
        for j = 1:3
            if( conv_flags(ext_loop*3 + j) == 0)
                if( mean_MSD_mat_comp(i, ext_loop*3 + j) <= SS_MSD_threshold )
                    conv_time_vec(ext_loop*3 + j) = i;
                    conv_flags(ext_loop*3 + j) = 1;
                end
            end
        end
    end
    
    subplot(1,3,ext_loop+1)
    plot(t,10*log10(mean_MSD_mat_comp(t,ext_loop*3 + 1)),t,10*log10(mean_MSD_mat_comp(t,ext_loop*3 + 2)),...
        t,10*log10(mean_MSD_mat_comp(t,ext_loop*3 + 3)), t, 10*log10(SS_MSD_threshold)*ones(size(t)),'r--' )
    
end

%% ===================================================================== %%
% ======================================================================= %
% == Generating Lines to be Used in LaTeX Table ========================= %
% ======================================================================= %

disp('======================================================')
inline_table = [];
for i = 1:9
    switch(i)
        case 1
            line_string = ['(L.I) &  LMS & ' num2str(alg_param_vec(i)) ' & (1) & ' ];
        case 2
            line_string = ['(R.I) &  RLS & ' num2str(alg_param_vec(i)) ' & (1) & ' ];
        case 3
            line_string = ['(N.I) &  LMS & ' num2str(alg_param_vec(i)) ' & (1) & ' ];
        case 4
            line_string = ['(L.II) &  LMS & ' num2str(alg_param_vec(i)) ' & (3) & ' ];
        case 5
            line_string = ['(R.II) &  RLS & ' num2str(alg_param_vec(i)) ' & (3) & ' ];
        case 6
            line_string = ['(N.II) &  LMS & ' num2str(alg_param_vec(i)) ' & (3) & ' ];
        case 7
            line_string = ['(L.III) &  LMS & ' num2str(alg_param_vec(i)) ' & (3) & ' ];
        case 8
            line_string = ['(R.III) &  RLS & ' num2str(alg_param_vec(i)) ' & (3) & ' ];
        case 9
            line_string = ['(N.III) &  LMS & ' num2str(alg_param_vec(i)) ' & (3) & ' ];
    end
    
    line_string = [ line_string num2str(SteadyState_MSD_vec(i)) ' & ' num2str(conv_time_vec(i)) ' & ' ...
                    num2str(mean( mean_elapsedTime_mat_comp(:,i))*10^(6)) ' \\ ' ];
                
   inline_table = [inline_table line_string ];
end
inline_table
disp('======================================================')
% ======================================================================= %

% ======================================================================= %
% ======================================================================= %