% ======================================================================= %
% == COPPE/UFRJ - Programa de Engenharia Eletrica (PEE) ================= %
% == Script: script01_Generate_Graph_Signal ============================= %
% == Responsible: Marcelo Jorge Mendes Spelta - Date: 2019/06/20 ======== %
% == E-mail: marcelo.spelta@smt.ufrj.br ================================= %
% ======================================================================= %

clc;clear;format shortEng;

% ======================================================================= %
% -- SETUP PARAMETERS FOR GRAPH GENERATION

K = 7;  % Number of considered stations next to the one desired (minimum number)
earthRadius = 6360; % km (approximation)

% ======================================================================= %
% -- Reading data from excel file --------------------------------------- %
[numbDataTable, textDataTable, rawDataTable] = ...
    xlsread("NOAA_Weather_Stations-Temperature_Jan2010",2);
% -------------------------------------
numberStations = size(numbDataTable,1);
% -------------------------------------

latitudeVector = numbDataTable(:,1); 
longitudeVector = numbDataTable(:,2); 

coordinates_matrix = [longitudeVector latitudeVector];   % 205x2-matrix containing the points longitude and latitude values

% ======================================================================= %

%% ===================================================================== %%
% ======================================================================= %
% == Generating adjacency matrix based on the Gaussian kernel and ======= %
% == the Haversine distance between nodes coordinates =================== %

% ======================================================================= %
% == DISTANCES MATRIX =================================================== %
% -- Haversine formula!
distancesMatrix = zeros( numberStations );
for extCounter = 1:(numberStations-1)
    for intCounter = (extCounter+1):numberStations
        % The next lines convert the rational degree values into radians for computing the Haversine distance
        x1 = longitudeVector(extCounter)*pi/180;
        y1 = latitudeVector(extCounter)*pi/180;
        x2 = longitudeVector(intCounter)*pi/180;
        y2 = latitudeVector(intCounter)*pi/180;
        
        distancesMatrix(extCounter, intCounter) = 2*earthRadius*asin( sqrt( sin((y2-y1)/2)^2 + cos(y1)*cos(y2)*sin((x2 - x1)/2)^2 ) );
        distancesMatrix(intCounter, extCounter) = distancesMatrix(extCounter, intCounter);
    end
end

A = zeros(numberStations);
meanDistance = mean(mean(distancesMatrix)); 
for extCounter = 1:(numberStations-1)
    for intCounter = (extCounter+1):numberStations
        A(extCounter, intCounter) = exp(-( (distancesMatrix(extCounter, intCounter)) / meanDistance));
        A(intCounter, extCounter) = A(extCounter, intCounter);
    end
end

% ======================================================================= %
% == INDICATION MATRIX ================================================== %

indicationMatrix = zeros( numberStations );

% -- Finding indicationMatrix based on the minimum distances between vertices
for stationCounter = 1:numberStations
    [ sorted_vec, sorted_vec_index] = sort( distancesMatrix( stationCounter, : ) );
    
    % The procedure included below is a very simple one for collecting at least $K$ neighbors for each graph node
    for intCounter = 2:(K+1) % I start taking values from 2 because it excludes the own node 
        indicationMatrix(stationCounter, sorted_vec_index(intCounter) ) = 1;
        indicationMatrix(sorted_vec_index(intCounter), stationCounter ) = 1;
    end
    
end

% ======================================================================= %
% == GENERATING ADJACENCY MATRIX
A = indicationMatrix.*A;

%% ===================================================================== %%
% == Obtaining temperature data for each weather station

% ======================================================================= %
% -- Reading data from excel file --------------------------------------- %
[numbDataTable2, textDataTable2, rawDataTable2] = ...
    xlsread("NOAA_Weather_Stations-Temperature_Jan2010",1);

day01_Hour01_Data = numbDataTable2(1:95:end,5);
day01_Hour12_Data = numbDataTable2(12:95:end,5);
day01_Hour18_Data = numbDataTable2(18:95:end,5);
save('Simu_Results/Illustrating_GSs','coordinates_matrix','A','day01_Hour01_Data','day01_Hour12_Data','day01_Hour18_Data')

completeDatasetMatrix = [];
for counter = 1:95
    completeDatasetMatrix = [completeDatasetMatrix numbDataTable2(counter:95:end,5)];
end

% For displaying variable: "plot(completeDatasetMatrix(100,:))"

%% ===================================================================== %%
% == For verifying if the temperature dataset can be considered a ======= %
% == bandlimited graph signal, we focus our attention to the particular = %
% == data from the month of July. ======================================= %

graph_signal =  day01_Hour01_Data;  % Temperatures at 01:00 January 1st - 2010     
% ======================================================================= %

L = diag(sum(A)) - A;     % Laplacian Matrix
[U,D] = eig(L);     	  % Eigendecomposition of the Adjacency matrix
freq = U'*graph_signal;   % Frequency-domain representation of the July data
[sorted_values, sorted_indices] = sort(abs(freq),'descend'); % Sorting frequency values in descending magnitude

% ======================================================================= %

for extCounter = 1:length(graph_signal)
    
    m = extCounter; % number of null components

    % -------------------------------------------
    used_indices = sort( sorted_indices(1:m), 'ascend' );
    U_f = [];
    s_f = [];
    
    for counter = 1:m
        U_f = [U_f U( : , used_indices(counter) ) ];
        s_f = [s_f; freq( used_indices(counter) )];
    end

    compact_graph_signal = U_f * s_f;
    % -----------------------------------------
    SD(extCounter) = (norm(compact_graph_signal - graph_signal, 2)) / (norm(graph_signal, 2)) ;
    % -----------------------------------------
    % -------------------------------------------
end

save('Simu_Results/SD_numberUsedComponents','SD')

% ======================================================================= %
% ======================================================================= %

%% ===================================================================== %%
% == For finally obtaining a bandlimited representation we consider ===== %
% == the F-largest frequency components of the GFT and construct the ==== %
% == matrix U_f. This matrix U_f, along with the number of sampled ====== %
% == nodes M, is used for obtaining the sampling set S and matrix D_s. == %
% ======================================================================= %
% == Both matrices U_f and D_s are stored since they will be used in ==== %
% == further simulations. =============================================== %
% ======================================================================= %

F = 125;        % number of used frequency components for the approximation

used_indices = sort( sorted_indices(1:F), 'ascend' );
U_f = [];
s_f = [];
    
for counter = 1:F
    U_f = [U_f U( : , used_indices(counter) ) ];
    s_f = [s_f; freq( used_indices(counter) )];
end

bandlimited_graph = U_f * s_f;  % Obtaining the approximated bandlimited graph signal based on U_f and s_f 

M = 130;

D_s = eig_sampling_strategy( M, U_f );

save('Simu_Results/General_Temperature_Data','D_s','U_f','completeDatasetMatrix');

% == END OF SCRIPT ====================================================== %
% ======================================================================= %