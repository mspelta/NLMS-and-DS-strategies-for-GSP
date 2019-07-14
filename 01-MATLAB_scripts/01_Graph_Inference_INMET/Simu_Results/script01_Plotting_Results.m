% ======================================================================= %
% == COPPE/UFRJ - Programa de Engenharia Eletrica (PEE) ================= %
% == Script: script01_Plotting_Results ================================== %
% == Responsible: Marcelo Jorge Mendes Spelta - Date: 2019/03/29 ======== %
% == E-mail: marcelo.spelta@smt.ufrj.br ================================= %
% ======================================================================= %

% ======================================================================= %
% == IMPORTANT: For plotting the graph pictures, it is necessary to run = %
% == the script "grasp_start" on folder "./GraSP-master/." ============== %

clc;clear;format shortEng;

addpath('../GraSP-master')   % Adding GraSP-master to MATLAB path
grasp_start                 % Initializing Grasp

%% ===================================================================== %%
% == Generating the empty Graph Signal figures ========================== %

clc

load('Illustrating_GSs')  % Loading graph signals related to the months of January, April and July

myGraph = grasp_plane_rnd( length(januaryData) );
myGraph.A = zeros(size(myGraph.A));
myGraph.layout = [coordinates_matrix(:,1), coordinates_matrix(:,2)];

myGraph.A_layout = 0;

x_offset = 2;
y_offset = 1;

limits = [ min(min( [januaryData aprilData julyData ] ) ) max(max( [januaryData aprilData julyData ] ) ) ];

%% ===================================================================== %%
% ----------------------------------------------------------------------- %
% -- JANUARY ------------------------------------------------------------ %
figure;
graph_signal =  januaryData; 
grasp_show_graph(gca, myGraph, 'node_values', graph_signal); 
xlim([min(coordinates_matrix(:,1))-x_offset max(coordinates_matrix(:,1))+x_offset ])
ylim([min(coordinates_matrix(:,2))-y_offset max(coordinates_matrix(:,2))+y_offset ])
xlabel('Longitude')
ylabel('Latitude')
colorbar
caxis(limits)

axis([min(coordinates_matrix(:,1))-x_offset max(coordinates_matrix(:,1))+x_offset min(coordinates_matrix(:,2))-y_offset max(coordinates_matrix(:,2))+y_offset ] )

% -- Generating January empty figure 
sizeFontImage = 40;
figProp = struct('size',sizeFontImage,'font','Times','lineWidth',2,'figDim',[1 1 800 800]);
figFileName = '../Figures/January_Null_Graph';
formatFig(gcf,figFileName,'en',figProp);

%% ===================================================================== %%
% ----------------------------------------------------------------------- %
% -- APRIL -------------------------------------------------------------- %
figure;
graph_signal =  aprilData; 
grasp_show_graph(gca, myGraph, 'node_values', graph_signal);
xlim([min(coordinates_matrix(:,1))-x_offset max(coordinates_matrix(:,1))+x_offset ])
ylim([min(coordinates_matrix(:,2))-y_offset max(coordinates_matrix(:,2))+y_offset ])
xlabel('Longitude')
ylabel('Latitude')
colorbar
caxis(limits)

axis([min(coordinates_matrix(:,1))-x_offset max(coordinates_matrix(:,1))+x_offset min(coordinates_matrix(:,2))-y_offset max(coordinates_matrix(:,2))+y_offset ] )

% -- Generating April empty figure 
figProp = struct('size',sizeFontImage,'font','Times','lineWidth',2,'figDim',[1 1 800 800]);
figFileName = '../Figures/April_Null_Graph';
formatFig(gcf,figFileName,'en',figProp);

%% ===================================================================== %%
% ----------------------------------------------------------------------- %
% -- JULY --------------------------------------------------------------- %
figure
graph_signal =  julyData; 
grasp_show_graph(gca, myGraph, 'node_values', graph_signal); 
xlim([min(coordinates_matrix(:,1))-x_offset max(coordinates_matrix(:,1))+x_offset ])
ylim([min(coordinates_matrix(:,2))-y_offset max(coordinates_matrix(:,2))+y_offset ])
xlabel('Longitude')
ylabel('Latitude')
colorbar
caxis(limits)

axis([min(coordinates_matrix(:,1))-x_offset max(coordinates_matrix(:,1))+x_offset min(coordinates_matrix(:,2))-y_offset max(coordinates_matrix(:,2))+y_offset ] )

% -- Generating July Empty figure 
figProp = struct('size',sizeFontImage,'font','Times','lineWidth',2,'figDim',[1 1 800 800]);
figFileName = '../Figures/July_Null_Graph';
formatFig(gcf,figFileName,'en',figProp);

%% ===================================================================== %%
% == Generating the complete Graph Signal figures for the months of ===== %
% == January, April and July. Complete in this case means the graph ===== %
% == signal along the graph structure. ================================== %
% ======================================================================= %

myGraph = grasp_plane_rnd( length(januaryData) );
myGraph.A = A;
myGraph.A_layout = A; 
myGraph.layout = [coordinates_matrix(:,1), coordinates_matrix(:,2)];

for extCounter = 1:length(januaryData)
    for intCounter = 1:length(januaryData)
        myGraph.distances(extCounter,intCounter) = 0;
    end
end

% ======================================================================= %
% ----------------------------------------------------------------------- %
% -- JANUARY ------------------------------------------------------------ %

figure
graph_signal =  januaryData;  
grasp_show_graph(gca, myGraph, 'node_values', graph_signal); 
xlim([min(coordinates_matrix(:,1))-x_offset max(coordinates_matrix(:,1))+x_offset ])
ylim([min(coordinates_matrix(:,2))-y_offset max(coordinates_matrix(:,2))+y_offset ])
xlabel('Longitude')
ylabel('Latitude')
colorbar
caxis(limits)

% -- Generating January complete figure 
figProp = struct('size',40,'font','Times','lineWidth',1,'figDim',[1 1 600 1000]);
figFileName = '../Figures/Graph_Signal_Months_January';
formatFig(gcf,figFileName,'en',figProp);

% ======================================================================= %
% ----------------------------------------------------------------------- %
% -- APRIL -------------------------------------------------------------- %

figure
graph_signal =  aprilData;
grasp_show_graph(gca, myGraph, 'node_values', graph_signal); 
xlim([min(coordinates_matrix(:,1))-x_offset max(coordinates_matrix(:,1))+x_offset ])
ylim([min(coordinates_matrix(:,2))-y_offset max(coordinates_matrix(:,2))+y_offset ])
colorbar
caxis(limits)

% -- Generating April complete figure 
figProp = struct('size',sizeFontImage,'font','Times','lineWidth',1,'figDim',[1 1 600 1000]);
figFileName = '../Figures/Graph_Signal_Months_April';
formatFig(gcf,figFileName,'en',figProp);

% ======================================================================= %
% ----------------------------------------------------------------------- %
% -- JULY --------------------------------------------------------------- %
figure
graph_signal =  julyData; 
grasp_show_graph(gca, myGraph, 'node_values', graph_signal);
xlim([min(coordinates_matrix(:,1))-x_offset max(coordinates_matrix(:,1))+x_offset ])
ylim([min(coordinates_matrix(:,2))-y_offset max(coordinates_matrix(:,2))+y_offset ])
colorbar
caxis(limits)

% -- Generating July complete figure
figProp = struct('size',sizeFontImage,'font','Times','lineWidth',1,'figDim',[1 1 600 1000]);
figFileName = '../Figures/Graph_Signal_Months_July';
formatFig(gcf,figFileName,'en',figProp);

%% ===================================================================== %%
% == Plotting the relation between reconstruction error and the number == %
% == of used frequency components. ====================================== %
% ======================================================================= %

load('SD_numberUsedComponents')

figure
k = 50:250;
plot(k, 100*SD(k), '-', 50:250 , 2.5*ones(201,1), 'r--' )
xlabel('Number of used components','Interpreter','latex','FontSize',18)
ylabel('Error ($\%$)','Interpreter','latex','FontSize',18)
grid on
axis([50 250 0 13] )
figProp = struct('size',20,'font','Times','lineWidth',4,'figDim',[1 1 800 350]); %1080]);%600 400]);
figFileName = '../Figures/Reconst_Error_Present';
formatFig(gcf,figFileName,'en',figProp);


%% ===================================================================== %%
% == Displaying the resemblance between the original July graph signal == %
% == and its bandlimited approximation. ================================= %

load('General_Bandlimited_GS_Data');

% ======================================================================= %
% ======================================================================= %

t = 1:length(graph_signal);
figure
plot(t,graph_signal,'b',t,bandlimited_graph,'r')
legend('Graph Signal','Bandlimited Graph Signal')

error_vector = (graph_signal - bandlimited_graph);
MSE = (error_vector' * error_vector ) / ( length(error_vector) );
MSD = norm(error_vector,2)/norm(graph_signal,2);

limits = [min(min(graph_signal), min(bandlimited_graph)) max(max(graph_signal), max(bandlimited_graph))];
x_offset = 2;
y_offset = 1;

figure
subplot(1,2,1)
grasp_show_graph(gca, myGraph, 'node_values', graph_signal); 
xlim([min(coordinates_matrix(:,1))-x_offset max(coordinates_matrix(:,1))+x_offset ])
ylim([min(coordinates_matrix(:,2))-y_offset max(coordinates_matrix(:,2))+y_offset ])
xlabel('Longitude')
ylabel('Latitude')
title_string = 'Graph Signal';
title(title_string)
colorbar
caxis(limits) 

subplot(1,2,2)
grasp_show_graph(gca, myGraph, 'node_values', bandlimited_graph); 
xlim([min(coordinates_matrix(:,1))-x_offset max(coordinates_matrix(:,1))+x_offset ])
ylim([min(coordinates_matrix(:,2))-y_offset max(coordinates_matrix(:,2))+y_offset ])
xlabel('Longitude')
ylabel('Latitude')
title_string = 'Bandlimited Graph Signal';
title(title_string)
colorbar
caxis(limits) 

% == END OF SCRIPT ====================================================== %
% ======================================================================= %