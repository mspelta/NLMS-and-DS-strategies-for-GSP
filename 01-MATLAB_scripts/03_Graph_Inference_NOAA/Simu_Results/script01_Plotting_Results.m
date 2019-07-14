% ======================================================================= %
% == COPPE/UFRJ - Programa de Engenharia Eletrica (PEE) ================= %
% == Script: script01_Plotting_Results ================================== %
% == Responsible: Marcelo Jorge Mendes Spelta - Date: 2019/06/20 ======== %
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

load('Illustrating_GSs')  % Loading graph signals related to 01:00, 12:00, and 18:00 January 1st 2010 

myGraph = grasp_plane_rnd( length(day01_Hour01_Data) );
myGraph.A = zeros(size(myGraph.A));
myGraph.layout = [coordinates_matrix(:,1), coordinates_matrix(:,2)];

myGraph.A_layout = 0;

x_offset = 2;
y_offset = 1;

limits = [ min(min( [day01_Hour01_Data day01_Hour12_Data day01_Hour18_Data ] ) ) ...
    max(max( [ day01_Hour01_Data day01_Hour12_Data day01_Hour18_Data ] ) ) ];

%% ===================================================================== %%
% ----------------------------------------------------------------------- %
% -- 01:00 January 1st -------------------------------------------------- %
figure;
graph_signal = day01_Hour01_Data; 
myGraph.A = A;
grasp_show_graph(gca, myGraph, 'node_values', graph_signal); 
xlim([min(coordinates_matrix(:,1))-x_offset max(coordinates_matrix(:,1))+x_offset ])
ylim([min(coordinates_matrix(:,2))-y_offset max(coordinates_matrix(:,2))+y_offset ])
xlabel('Longitude')
ylabel('Latitude')
colorbar
caxis(limits)

axis([min(coordinates_matrix(:,1))-x_offset max(coordinates_matrix(:,1))+x_offset ...
    min(coordinates_matrix(:,2))-y_offset max(coordinates_matrix(:,2))+y_offset ] )

% -- Generating 01:00 January 1st's empty figure 
sizeFontImage = 40;
figProp = struct('size',sizeFontImage,'font','Times','lineWidth',2,'figDim',[1 1 800 500]);
figFileName = '../Figures/day01_Hour01_Data';
formatFig(gcf,figFileName,'en',figProp);

%% ===================================================================== %%
% ----------------------------------------------------------------------- %
% -- 12:00 January 1st -------------------------------------------------- %
figure;
graph_signal = day01_Hour12_Data; 
grasp_show_graph(gca, myGraph, 'node_values', graph_signal);
xlim([min(coordinates_matrix(:,1))-x_offset max(coordinates_matrix(:,1))+x_offset ])
ylim([min(coordinates_matrix(:,2))-y_offset max(coordinates_matrix(:,2))+y_offset ])
xlabel('Longitude')
ylabel('Latitude')
colorbar
caxis(limits)

axis([min(coordinates_matrix(:,1))-x_offset max(coordinates_matrix(:,1))+x_offset ...
    min(coordinates_matrix(:,2))-y_offset max(coordinates_matrix(:,2))+y_offset ] )

% -- Generating 12:00 January 1st's empty figure 
figProp = struct('size',sizeFontImage,'font','Times','lineWidth',2,'figDim',[1 1 800 500]);
figFileName = '../Figures/day01_Hour12_Data';
formatFig(gcf,figFileName,'en',figProp);

%% ===================================================================== %%
% ----------------------------------------------------------------------- %
% -- 18:00 January 1st -------------------------------------------------- %
figure
graph_signal = day01_Hour18_Data;
grasp_show_graph(gca, myGraph, 'node_values', graph_signal); 
xlim([min(coordinates_matrix(:,1))-x_offset max(coordinates_matrix(:,1))+x_offset ])
ylim([min(coordinates_matrix(:,2))-y_offset max(coordinates_matrix(:,2))+y_offset ])
xlabel('Longitude')
ylabel('Latitude')
colorbar
caxis(limits)

axis([min(coordinates_matrix(:,1))-x_offset max(coordinates_matrix(:,1))+x_offset ...
    min(coordinates_matrix(:,2))-y_offset max(coordinates_matrix(:,2))+y_offset ] )

% -- Generating 18:00 January 1st's empty figure  
figProp = struct('size',sizeFontImage,'font','Times','lineWidth',2,'figDim',[1 1 800 500]);
figFileName = '../Figures/day01_Hour18_Data';
formatFig(gcf,figFileName,'en',figProp);

% == END OF SCRIPT ====================================================== %
% ======================================================================= %