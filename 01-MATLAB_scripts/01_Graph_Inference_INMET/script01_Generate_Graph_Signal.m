% ======================================================================= %
% == COPPE/UFRJ - Programa de Engenharia Eletrica (PEE) ================= %
% == Script: script01_Generate_Graph_Signal ============================= %
% == Responsible: Marcelo Jorge Mendes Spelta - Date: 2019/03/29 ======== %
% == E-mail: marcelo.spelta@smt.ufrj.br ================================= %
% ======================================================================= %

clc;clear;format shortEng;

% ======================================================================= %
% -- SETUP PARAMETERS FOR GRAPH GENERATION

K = 8;  % Number of considered stations next to the one desired (minimum number)
earthRadius = 6360; % km (approximation)

% ======================================================================= %
% -- Reading data from excel file --------------------------------------- %
[numbDataTable, textDataTable, rawDataTable] = ...
    xlsread("Brazilian_Weather_Stations-Temperature_1961-1990",1);
% -------------------------------------
numberStations = size(numbDataTable,1);
% -------------------------------------
% It adjust textDataTable because the first line is presented!
textDataTable(1:numberStations, :) = textDataTable(2:(numberStations+1), :);
textDataTable(numberStations+1, :) = [];
% -------------------------------------

myDataTable = zeros( numberStations, 17);
for vert_counter = 1:numberStations
    myDataTable(vert_counter, 1) = numbDataTable(vert_counter, 1);
    
    myDataTable(vert_counter, 4:16) = numbDataTable(vert_counter, 7:19 ); 
end

for hor_counter = 1:1:2
    for vert_counter = 1:numberStations
        posRepresentation = textDataTable(vert_counter, 3 + hor_counter);
        charPosRepres = char( posRepresentation );
        % The following command turns the Degree Minute representation in a rational degree representation
        myDataTable(vert_counter, 1 + hor_counter) = degMin_2_deg(str2num(charPosRepres(1:2)),str2num(charPosRepres(4:5))); % dms2degrees( [str2num(charPosRepres(1:2)) str2num(charPosRepres(4:5)) 0] );
        
        if ( ( charPosRepres(7) == 'S' ) || ( charPosRepres(7) == 'W' ) ) % Converting south and west directions into negative numbers
            myDataTable(vert_counter, 1 + hor_counter) = - myDataTable(vert_counter, 1 + hor_counter);
        end
            
    end
end

latitudeVector = myDataTable(:, 2); 
longitudeVector = myDataTable(:, 3); 

coordinates_matrix = [longitudeVector latitudeVector];   % 299x2-matrix containing the points longitude and latitude values

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

januaryData = myDataTable(:,4);
aprilData = myDataTable(:,7);
julyData = myDataTable(:,10);
save('Simu_Results/Illustrating_GSs','coordinates_matrix','A','januaryData','aprilData','julyData')

%% ===================================================================== %%
% == For verifying if the temperature dataset can be considered a ======= %
% == bandlimited graph signal, we focus our attention to the particular = %
% == data from the month of July. ======================================= %

graph_signal =  myDataTable(:, 10);     % Data related to the month of July    
% ======================================================================= %

[U,D] = eig(A);     	  % Eigendecomposition of the Adjacency matrix
freq = U'*graph_signal;   % Frequency-domain representation of the July data
[sorted_values, sorted_indices] = sort(abs(freq),'descend'); % Sorting frequency values in descending magnitude

% ======================================================================= %

for extCounter = 1:length(graph_signal)
    
    m = extCounter; %extCounter*75; %150;    % number of null components

    % -------------------------------------------
    %indication_vector = zeros(length(sorted_values), 1);
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

F = 200; % number of used frequency components for the approximation

used_indices = sort( sorted_indices(1:F), 'ascend' );
U_f = [];
s_f = [];
    
for counter = 1:F
    U_f = [U_f U( : , used_indices(counter) ) ];
    s_f = [s_f; freq( used_indices(counter) )];
end

bandlimited_graph = U_f * s_f;  % Obtaining the approximated bandlimited graph signal based on U_f and s_f 

M = 210;

D_s = eig_sampling_strategy( M, U_f );

save('Simu_Results/General_Bandlimited_GS_Data','D_s','U_f','graph_signal','bandlimited_graph');

% == END OF SCRIPT ====================================================== %
% ======================================================================= %