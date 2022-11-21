%%%%%%%%%% THIS IS FINAL PHELPS MODEL CODE

clc,clear all;

%% READING THE PROBE LIST POINTS

ProbeList = readmatrix('D:\TU Delft\TNO Thesis\Final Report\Final Probes X section\ProbeList_InsertC117.csv');
ProbeList(1,:) = []; % TO DELETE THE FIRST ROW OF THE READ FILE
ProbeList(:,1) = []; % TO DELETE THE FIRST COLUMN OF THE READ FILE

% Read simulations data
ResultDataTable = readtable('D:\TU Delft\TNO Thesis\Final Report\Final Probes X section\ResultData_InsertC117.csv');


%%%%%%%%%%%%%%%%% PRE-PROCESSING OF THE SIMULATION DATA STARTS HERE %%%%%%%%%%%%%%%%%%%%%%%%%% 
%% NOW WE PROCESS ALL THE DATA EXTRACTED FROM MOLDEX3D
% We have the data as TABLE for each of the point.

probes = length(ProbeList(:,1)); % No. Of probes:

% Stores all the times for the nodes
CellTime = cell(1,probes);  
for i=1:probes
    column = ResultDataTable(:,2*i);
%     column = ResultDataMatrix(:,2*i);
    
    CellTime{1,i} = column; 
end



CellAverageFiberLengthLW = cell(1,probes);  
indexAF = 3;   % This is the first index of fiber length by number. This needs to be manually checked.                             

for i=1:probes
    column = ResultDataTable(:,indexAF);  
    CellAverageFiberLengthLW{1,i} = column; 
    indexAF = indexAF + 2;
end


CellAverageFiberLengthLN = cell(1,probes); 
indexAF2 = 3+(2*probes);                               


for i=1:probes
    column = ResultDataTable(:,indexAF2);  
    CellAverageFiberLengthLN{1,i} = column; 
    indexAF2 = indexAF2 + 2;
end



CellShearRate = cell(1,probes); 
indexSR = 3+(2*2*probes);        


for i=1:probes
    column = ResultDataTable(:,indexSR); 
    CellShearRate{1,i} = column;
    indexSR = indexSR + 2;
end


CellViscosity = cell(1,probes);  
indexV = 3+(3*2*probes);       

for i=1:probes
    column = ResultDataTable(:,indexV);  
    CellViscosity{1,i} = column;
    indexV = indexV + 2;
end

CellVelocityX = cell(1,probes); 
indexVelX = 3+(4*2*probes);      

for i=1:probes
    column = ResultDataTable(:,indexVelX);  
    CellVelocityX{1,i} = column;
    indexVelX = indexVelX + 2;
end


CellVelocityY = cell(1,probes); 
indexVelY = 3+(5*2*probes);   

for i=1:probes
    column = ResultDataTable(:,indexVelY);
    CellVelocityY{1,i} = column;
    indexVelY = indexVelY + 2;
end


CellVelocityZ = cell(1,probes); 
indexVelZ = 3+(6*2*probes);      

for i=1:probes
    column = ResultDataTable(:,indexVelZ);   
    CellVelocityZ{1,i} = column;
    indexVelZ = indexVelZ + 2;
end

CellVelocityTot = cell(1,probes);  
indexVelTot = 3+(7*2*probes);      

for i=1:probes
    column = ResultDataTable(:,indexVelTot); 
    CellVelocityTot{1,i} = column;
    indexVelTot = indexVelTot + 2;
end


%% Now we link the probes with Cells with result data that we have created. 

% Maybe create a cell to save the coordinates
CellCoordinate = cell(1,3);  % If we want to specify size

for i=1:3
    column = ProbeList(:,i); 
    CellCoordinate{1,i} = column;
end

%% We find the coordinate to start the simulation with.
% we can find the index of the CellTime with largest length. We can do this
% by checking the values of each column in the CellTime and look for NaN
% values. Then the cell column with least Nan values is has most times. And
% Corresponding point is the starting point. 

counter = zeros(1,probes); % This list stores the no. of Nan for each of the probes points 

for i = 1:length(CellTime) 
    %Check for Nan and add counter    
    for j=1:height(CellTime{1,i})   % if TABLE Cell is used

        celltimearray = table2array(CellTime{1,i});
        if isnan(celltimearray(j)) 
            counter(i) = counter(i) + 1;    
        end
    end 
    
    %%% THIS PART OF THE CODE IS TO ORGANISE THE DATA SUCH THAT NaN VALUES
    %%% ARE ON TOP BEFORE WE GET PROPER READINGS.
    CellTime{1,i} = circshift(CellTime{1,i},counter(i));
    CellAverageFiberLengthLW{1,i} = circshift(CellAverageFiberLengthLW{1,i},counter(i));
    CellAverageFiberLengthLN{1,i} = circshift(CellAverageFiberLengthLN{1,i},counter(i));
    CellShearRate{1,i} = circshift(CellShearRate{1,i},counter(i));
    CellViscosity{1,i} = circshift(CellViscosity{1,i},counter(i));
    CellVelocityTot{1,i} = circshift(CellVelocityTot{1,i},counter(i));
    CellVelocityX{1,i} = circshift(CellVelocityX{1,i},counter(i));
    CellVelocityY{1,i} = circshift(CellVelocityY{1,i},counter(i));
    CellVelocityZ{1,i} = circshift(CellVelocityZ{1,i},counter(i));
end

% NOW USING THE COUNTER, WE SHIFT THE TABLE COLUMNS FROM ALL THE RESULTS. 

% Create an array to store the probe IDs
nodeID = zeros(1,probes); % We put probes-1, because we dont need probe 1
for i=1:probes
    nodeID(i) = i;
end


[counterNew,sortidx] = sort(counter,'ascend'); % To sort the probes times in ascending order of Nan values. 
nodeIDSort = nodeID(sortidx); % We have here created the order of the nodes in which they get filled up.

%%% MAYBE ITS BETTER NOW TO SORT ALL THE CELL DATA IN THIS ORDER AND MOVE
%%% AHEAD. 
CellTimeSort = CellTime(sortidx);
CellAverageFiberLengthLWSort = CellAverageFiberLengthLW(sortidx); 
CellAverageFiberLengthLNSort = CellAverageFiberLengthLN(sortidx);
CellShearRateSort = CellShearRate(sortidx);
CellViscositySort = CellViscosity(sortidx);
CellVelocityXSort = CellVelocityX(sortidx);
CellVelocityYSort = CellVelocityY(sortidx);
CellVelocityZSort = CellVelocityZ(sortidx);
CellVelocityTotSort = CellVelocityTot(sortidx);
CellCoordinateSort = cell(1,3);  % If we want to specify size
for i=1:3
    CellCoordinateSort{1,i} = CellCoordinate{1,i}(sortidx);  
end

%%% WE SHOULD HERE ALSO DELETE THE NaN VALUES FROM THE FIRST FILLED
%%% POINT AND THE CORRESPONDING TABLES

if counterNew(1) > 0
    for i = 1:length(CellTime) 
        % TODO
        CellTimeSort{1,i}(1:counterNew(1),:) = [];
        CellAverageFiberLengthLWSort{1,i}(1:counterNew(1),:) = [];
        CellAverageFiberLengthLNSort{1,i}(1:counterNew(1),:) = [];
        CellShearRateSort{1,i}(1:counterNew(1),:) = [];
        CellViscositySort{1,i}(1:counterNew(1),:) = [];
        CellVelocityXSort{1,i}(1:counterNew(1),:) = [];
        CellVelocityYSort{1,i}(1:counterNew(1),:) = [];
        CellVelocityZSort{1,i}(1:counterNew(1),:) = [];
        CellVelocityTotSort{1,i}(1:counterNew(1),:) = [];
    end
end

%% IN THIS SECTION WE TRY TO CREATE STRUCTURE TO STORE ALL THE VALUES

CellNodeID = cell(1,length(nodeIDSort));
CellX = cell(1,length(nodeIDSort));     
CellY = cell(1,length(nodeIDSort));
CellZ = cell(1,length(nodeIDSort));

for i=1:length(nodeIDSort) 
    CellNodeID(1,i) = {nodeIDSort(i)};
    CellX(1,i) = {CellCoordinateSort{1,1}(i)};
    CellY(1,i) = {CellCoordinateSort{1,2}(i)};
    CellZ(1,i) = {CellCoordinateSort{1,3}(i)};
end

RunData = struct('NodeID',CellNodeID,'Time',CellTimeSort,'FiberLengthLW',CellAverageFiberLengthLWSort,'FiberLengthLN',CellAverageFiberLengthLNSort,'ShearRate',CellShearRateSort,...
    'Viscosity',CellViscositySort,'VelocityX',CellVelocityXSort,'VelocityY',CellVelocityYSort,'VelocityZ',CellVelocityZSort,'XPoint',CellX,'YPoint',CellY,'ZPoint',CellZ);


%%  Lets define the grid
FieldPoints = sum(arrayfun(@(RunData) ~isempty(RunData.XPoint),RunData)) ; % Calculates the number of points 
                                                                           % in the field. Should be equal to variable 'probes'
for i = 1:FieldPoints
    Xtemp(i) = RunData(i).XPoint;
    Ytemp(i) = RunData(i).YPoint;
    Ztemp(i) = RunData(i).ZPoint;
end
                                                                           

X1 = length(unique(Xtemp));%No. of Grid points in X direction
Y1 = length(unique(Ytemp));%No. of Grid points in Y direction
Z1 = length(unique(Ztemp));%No. of Grid points in Z direction

%%%% HERE WE NEED TO MAKE SURE THAT THE GENERATED DIRECTION OF THE x1,x2,x3
%%%% IS IN THE FLOW DIRECTION. THEREFORE, WE MAY OR MAY NOT NEED TO FLIP
%%%% THE RESULTING ARRAY. 

x1 = unique(Xtemp); % distance point in x direction
y1 = unique(Ytemp); % distance point in y direction
z1 = unique(Ztemp); % distance point in z direction

x1 = flip(x1); % Optional --- if the flow direction is in negative direction of the axis
% y1 = flip(y1);
% z1 = flip(z1);

%%  Now lets add coordinates to the structure:
CellXcoord = cell(1,length(nodeIDSort));    % To store the x cordinate in from discretization
CellYcoord = cell(1,length(nodeIDSort));    % To store the y cordinate in from discretization
CellZcoord = cell(1,length(nodeIDSort));    % To store the z cordinate in from discretization

for i = 1:FieldPoints
    CellXcoord(1,i) = {find(x1 == RunData(i).XPoint)};
    CellYcoord(1,i) = {find(y1 == RunData(i).YPoint)};
    CellZcoord(1,i) = {find(z1 == RunData(i).ZPoint)};
end

%%% Structure to store all the data in an organised manner.
RunData = struct('NodeID',CellNodeID,'Time',CellTimeSort,'FiberLengthLW',CellAverageFiberLengthLWSort,'FiberLengthLN',CellAverageFiberLengthLNSort,...
    'ShearRate',CellShearRateSort,...
    'Viscosity',CellViscositySort,'VelocityX',CellVelocityXSort,'VelocityY',CellVelocityYSort,'VelocityZ',CellVelocityZSort,'XPoint',CellX,'YPoint',CellY,'ZPoint',CellZ,...
    'XCoord',CellXcoord,'YCoord',CellYcoord,'ZCoord',CellZcoord);

%% We write the code here to find the layers in each of the directions of the geometry

layerX =  length(unique(Xtemp)); % Obtained due to limited sections in X direction.
sample1 = x1(1); % Picking the first element as sample in X 

repx = 0;
while sample1 == Xtemp(repx+1)
    repx = repx+1;
end

% For finding layers in y:
layerY = 0;

for repy = 1:repx
    if (Xtemp(repy) == Xtemp(1)) && (Ztemp(repy) == Ztemp(1))
        layerY = layerY+1;
    end
end

% For finding layers in z:
layerZ = 0;

for repz = 1:repx
    if (Xtemp(repz) == Xtemp(1)) && (Ytemp(repz) == Ytemp(1))
        layerZ = layerZ+1;
    end
end

layer = repx; % as the number of total nodes in each x section is equal to repx. 

%%%%%%%%%%%%%%%%% PRE-PROCESSING OF THE SIMULATION DATA ENDS HERE %%%%%%%%%%%%%%%%%%%%%%%%%% 

%%%%%%%%%%%%%%%%% MODEL APPLICATION OF THE SIMULATION DATA STARTS HERE %%%%%%%%%%%%%%%%%%%%%%%%%% 

%% Now we go ahead with rest of the logic for simulation with Bechara Model.

%%%% ASSIGNING MATERIALS PROPERTIES FOR USE 

% Initiate the constants:
Ef = 72E9; % Pa Young's modulus of the fiber 
sig_ut = 2.6E9; % Pa Ultimate tensile strength of the fiber 
df = 12*10^-6; % m fiber diameter

%%% PLACING THE PARAMETERS HERE FROM THE BECHARA PAPER %%%%
lambda = 0.085; % For 30% fiber concentration (taken from bechara paper)
xi1Lw = 12.1E-4;
xi2Lw = 4.66E-4; % We should use pre-dispersed coefficients because there is breakage in the screws. 

xi1Ln = 14.5E-4;
xi2Ln = 4.97E-4;

lmax = 3*10^-3 ; % m fiber length as input

timeFilling = [0; table2array(CellTimeSort{1,1})]; % Defining the time array
tstep = length(timeFilling);        % Defining time step

%%%%%%% FILLING THE INITIAL POINTS
temp = size(CellTimeSort);
ProbePoints = temp(2);
save_index_first = []; % THE FIRST INDICATES THAT THIS IS THE INITIAL POINT WHERE FIBERS ARE INTRODUCED.

for i=1:ProbePoints
    if ismember(timeFilling(2),table2array(CellTimeSort{1,i})) % WE PUT THE INDEX OF timeFilling AS 2 IF FIRST VAL IS 0.
        save_index_first = [save_index_first, i];     
    end
end

% Now we fill the corresponding coordinates in the N matrix. To do that,
% get the point from the index in save_index. 


for i = 1:length(save_index_first)
    xpt(i) = RunData(save_index_first(i)).XPoint;
    ypt(i) = RunData(save_index_first(i)).YPoint;
    zpt(i) = RunData(save_index_first(i)).ZPoint;
end

for i = 1:length(save_index_first)
    idx1_first(i) = RunData(save_index_first(i)).XCoord;
    idx2_first(i) = RunData(save_index_first(i)).YCoord;
    idx3_first(i) = RunData(save_index_first(i)).ZCoord;
end


% Now assining the average fiber lengths (Lw/Ln) at the first probe equal to input. 
LN = zeros(1,tstep,X1,Y1,Z1); % FOR 3D case
LW = zeros(1,tstep,X1,Y1,Z1); % FOR 3D case

for i = 1:length(idx1_first)

    LN(1,1,idx1_first(i),idx2_first(i),idx3_first(i)) = lmax; % Assigning the initial points with complete fiber lengths 
    LW(1,1,idx1_first(i),idx2_first(i),idx3_first(i)) = lmax; % Assigning the initial points with complete fiber lengths
    
    % OR Providing variable fiber length at the input:
%     LN(1,1,idx1_first(i),idx2_first(i),idx3_first(i)) = RunData(i).FiberLengthLW{1,1}*1E-3; % Assigning the initial points with varying fiber lengths 
%     LW(1,1,idx1_first(i),idx2_first(i),idx3_first(i)) = RunData(i).FiberLengthLN{1,1}*1E-3; % Assigning the initial points with varying fiber lengths
   
                               
end

for T = 2:tstep

    delta_t = timeFilling(T)-timeFilling(T-1); 
    
    % Check again if the probe(s) is/are filled at the given time.
    save_index = [];
    idx1 = [];
    idx2 = [];
    idx3 = [];
    
    for i=1:ProbePoints
        if ismember(timeFilling(T),table2array(CellTimeSort{1,i})) % WE PUT THE INDEX OF timeFilling as T.
            save_index = [save_index, i];     
        end
    end

    for i = 1:length(save_index)
        xpt(i) = RunData(save_index(i)).XPoint;
        ypt(i) = RunData(save_index(i)).YPoint;
        zpt(i) = RunData(save_index(i)).ZPoint; 
    end

    for i = 1:length(save_index)
        idx1(i) = RunData(save_index(i)).XCoord;
        idx2(i) = RunData(save_index(i)).YCoord;
        idx3(i) = RunData(save_index(i)).ZCoord;
    end
    
    for m = 1:length(idx1)  

        %%%% For timesteps with initial as 0.               
        gamma = RunData(m).ShearRate{T-1,1};        
        eta_m = RunData(m).Viscosity{T-1,1}*0.1; % 0.1 is multiplied to change the units to SI units. 
        v1 = RunData(m).VelocityX{T-1,1}*0.01; % Velocity in the x direction, and change units to m/s
        v2 = RunData(m).VelocityY{T-1,1}*0.01; % Velocity in the y direction, and change units to m/s
        v3 = RunData(m).VelocityZ{T-1,1}*0.01; % Velocity in the z direction, and change units to m/s
        
%         %%%% For timesteps without initial as 0.
%         gamma = RunData(m).ShearRate{T,1};        
%         eta_m = RunData(m).Viscosity{T,1}*0.1; % 0.1 is multiplied to change the units to SI units. 
%         v1 = RunData(m).VelocityX{T,1}*0.01; % Velocity in the x direction, and change units to m/s 
%         v2 = RunData(m).VelocityY{T,1}*0.01; % Velocity in the y direction, and change units to m/s
%         v3 = RunData(m).VelocityZ{T,1}*0.01; % Velocity in the z direction, and change units to m/s


        % Calculate the Linf
        Linf = lambda*((sig_ut*df^2)/(eta_m*gamma)).^(1/2); 
        
        % Calculating breakage rate constants: 
        k1Lw = xi1Lw*gamma;
        k2Lw = xi2Lw*gamma;

        k1Ln = xi1Ln*gamma;
        k2Ln = xi2Ln*gamma;
        
        
        
        % INLET FACE FILLED AT BEGINNING
        if idx1(m) == 1 && (isnan(RunData(m).Time{1,1})) == 0
            LW(1,T,idx1(m),idx2(m),idx3(m)) = LW(1,T-1,idx1(m),idx2(m),idx3(m)) - k2Lw.*(LW(1,T-1,idx1(m),idx2(m),idx3(m)) - Linf)*delta_t;
            LN(1,T,idx1(m),idx2(m),idx3(m)) = LN(1,T-1,idx1(m),idx2(m),idx3(m)) - k2Ln.*(LN(1,T-1,idx1(m),idx2(m),idx3(m)) - Linf)*delta_t;

        % INLET FACE FILLED IN T2
        elseif idx1(m) == 1 && (isnan(RunData(m).Time{1,1}))
            LW(1,T,idx1(m),idx2(m),idx3(m)) = LW(1,T-1,idx1_first(1),idx2_first(1),idx3_first(1))- k2Lw.*(LW(1,T-1,idx1_first(1),idx2_first(1),idx3_first(1)) - Linf)*delta_t;
            LN(1,T,idx1(m),idx2(m),idx3(m)) = LN(1,T-1,idx1_first(1),idx2_first(1),idx3_first(1))- k2Ln.*(LN(1,T-1,idx1_first(1),idx2_first(1),idx3_first(1)) - Linf)*delta_t;
        elseif idx1(m)>1
                                   
            %%%% Connecting the correct nodes
            node = RunData(m).NodeID;
            node1 = node-layer;

            for ii = 1:FieldPoints
                if RunData(ii).NodeID == node1 && ii>m
                    node1 = node1-layer;
                    for jj = 1:FieldPoints
                        if RunData(jj).NodeID == node1
                            m1 = jj;
                            break
                        end
                    end
                elseif RunData(ii).NodeID == node1 && ii<m
                    m1 = ii;
                    break                     
                end
            end             

            % We calculate the distance. 
            distance = sqrt((RunData(m).XPoint*1E-3 - RunData(m1).XPoint*1E-3 ).^2 + (RunData(m).YPoint*1E-3 - RunData(m1).YPoint*1E-3 ).^2 +(RunData(m).ZPoint*1E-3 - RunData(m1).ZPoint*1E-3 ).^2 );
            xdistance = sqrt((RunData(m).XPoint*1E-3 - RunData(m1).XPoint*1E-3 ).^2 );
            ydistance = sqrt((RunData(m).YPoint*1E-3 - RunData(m1).YPoint*1E-3 ).^2 );
            zdistance = sqrt((RunData(m).ZPoint*1E-3 - RunData(m1).ZPoint*1E-3 ).^2);

 
                


%%%%%% UPDATING WITH THE SAME NODE:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if LW(1,T-1,idx1(m),idx2(m),idx3(m)) == 0  
                if LW(1,T-1,idx1(m1),idx2(m1),idx3(m1)) == 0
                    %TODO
                    LW(1,T,idx1(m),idx2(m),idx3(m)) = LW(1,T,idx1(m1),idx2(m1),idx3(m1));
                else
                    LW(1,T,idx1(m),idx2(m),idx3(m)) = LW(1,T-1,idx1(m1),idx2(m1),idx3(m1));
                end
            else
                % TODO
                LW(1,T,idx1(m),idx2(m),idx3(m)) = LW(1,T-1,idx1(m),idx2(m),idx3(m)) - k2Lw.*(LW(1,T-1,idx1(m1),idx2(m1),idx3(m1)) - Linf)*delta_t;
            end
            if LN(1,T-1,idx1(m),idx2(m),idx3(m)) == 0    
                if LN(1,T-1,idx1(m1),idx2(m1),idx3(m1)) == 0
                    %TODO
                    LN(1,T,idx1(m),idx2(m),idx3(m)) = LN(1,T,idx1(m1),idx2(m1),idx3(m1));
                else
                    LN(1,T,idx1(m),idx2(m),idx3(m)) = LN(1,T-1,idx1(m1),idx2(m1),idx3(m1));
                end
            else
                % TODO
                LN(1,T,idx1(m),idx2(m),idx3(m)) = LN(1,T-1,idx1(m),idx2(m),idx3(m)) - k2Ln.*(LN(1,T-1,idx1(m1),idx2(m1),idx3(m1)) - Linf)*delta_t;
            end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ENDS HERE !!!!!!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%         % Fitting to the modulus, taking away negative conditions:
%         LW(1,T,idx1(m),idx2(m),idx3(m)) = abs(LW(1,T,idx1(m),idx2(m),idx3(m)));
%         LN(1,T,idx1(m),idx2(m),idx3(m)) = abs(LN(1,T,idx1(m),idx2(m),idx3(m)));

        end          
        
    end
end

%% Calculating the fiber lengths at the corresponding points for plotting.

GraphArrayLW = zeros(FieldPoints,2); % Col 1: Moldex3d Col 2: Bechara 
GraphArrayLN = zeros(FieldPoints,2); % Col 1: Moldex3d Col 2: Bechara 
GraphArraySR = zeros(FieldPoints,1); % For storing shear rate
GraphArrayVis = zeros(FieldPoints,1); % For storing viscosities
% GraphArrayN = zeros(FieldPoints,M); % For storing FLD at end for each point

for i = 1:FieldPoints
    j = RunData(i).NodeID;
    
    
    
    GraphArrayLW(j,1) = RunData(i).FiberLengthLW{end,1}; 
    GraphArrayLN(j,1) = RunData(i).FiberLengthLN{end,1};
    GraphArraySR(j,1) = RunData(i).ShearRate{end,1};
    GraphArrayVis(j,1) = RunData(i).Viscosity{end,1}.*0.1;
    
    xx = RunData(i).XCoord;
    yy = RunData(i).YCoord;
    zz = RunData(i).ZCoord;
    
    %     if LW(1,end,xx,yy,zz) >= lmax
%         LW(1,end,xx,yy,zz) = 0.9*lmax;
%     end
%     
%     if LN(1,end,xx,yy,zz) >= lmax
%         LN(1,end,xx,yy,zz) = 0.9*lmax;
%     end
%     
%     if LW(1,end,xx,yy,zz) <= 0
%         LW(1,end,xx,yy,zz) = delta_l;
%     end
%     
%     if LN(1,end,xx,yy,zz) <= 0
%         LN(1,end,xx,yy,zz) = delta_l;
%     end
    
    
    GraphArrayLW(j,2) = LW(1,end,xx,yy,zz);

    
    GraphArrayLN(j,2) = LN(1,end,xx,yy,zz);

    
%     for mm = 1:M
%         GraphArrayN(j,mm) = Nstore(end,mm,xx,yy,zz);
%     end
    
    
end

%%% Smoothing the noise data:
SmoothArrayLWBechara = smoothdata(GraphArrayLW(:,2),'rloess');
SmoothArrayLNBechara = smoothdata(GraphArrayLN(:,2),'rloess');

figure()
plot(GraphArrayLW(:,1))
hold on 
plot(GraphArrayLW(:,2)*1000)
hold on 
plot(SmoothArrayLWBechara*1000,'g')
hold off
legend('Simulations Moldex3D','Bechara Model MATLAB','Smooth Bechara Model MATLAB')
xlabel('probes')
ylabel('fiber lengths (mm)')
title('measurements for each probe')


% % For section wise division: 
sections = ceil(probes/layer);
GraphSectionLW = zeros(sections,2); % Col 1: Moldex3d Col 2: Phelps 
GraphSectionLN = zeros(sections,2); % Col 1: Moldex3d Col 2: Phelps
SmoothSection = zeros(sections,2);% Col 1: PhelpsLW Col 2: PhelpsLN 

for kk = 1:sections
    
    if kk ~= sections
        GraphSectionLW(kk,1) = mean(GraphArrayLW((1+(kk-1)*layer):kk*layer,1));
        GraphSectionLW(kk,2) = mean(GraphArrayLW((1+(kk-1)*layer):kk*layer,2));

        GraphSectionLN(kk,1) = mean(GraphArrayLN((1+(kk-1)*layer):kk*layer,1));
        GraphSectionLN(kk,2) = mean(GraphArrayLN((1+(kk-1)*layer):kk*layer,2));

        SmoothSection(kk,1) = mean(SmoothArrayLWBechara((1+(kk-1)*layer):kk*layer,1));
        SmoothSection(kk,2) = mean(SmoothArrayLNBechara((1+(kk-1)*layer):kk*layer,1));
    else
        % TODO
        GraphSectionLW(kk,1) = mean(GraphArrayLW((1+(kk-1)*layer):end,1));
        GraphSectionLW(kk,2) = mean(GraphArrayLW((1+(kk-1)*layer):end,2));

        GraphSectionLN(kk,1) = mean(GraphArrayLN((1+(kk-1)*layer):end,1));
        GraphSectionLN(kk,2) = mean(GraphArrayLN((1+(kk-1)*layer):end,2));

        SmoothSection(kk,1) = mean(SmoothArrayLWBechara((1+(kk-1)*layer):end,1));
        SmoothSection(kk,2) = mean(SmoothArrayLNBechara((1+(kk-1)*layer):end,1));
    
    end


end

figure()
plot(GraphSectionLW(:,1))
hold on 
plot(GraphSectionLW(:,2)*1000)
title("\xi = " + xi2Lw +" \lambda = "+lambda)
xlabel('sections')
ylabel('Fiber Lengths (mm)')
legend('Simulations Moldex3D','Bechara Model MATLAB')
% xlim([2 7])
% xticks([3 4 5 6])



%% Main sections to be defined for separating fiber lengths in the parts. 

% Sections are x3 = -15.5 x4 = -25.5 x5 = -45.5
% col 1: part 3 col 2: part 4 col 3: part 5 col 4: part 6 
SectionsLWBechara = zeros(1,4); 
SectionsLNBechara = zeros(1,4);
% SectionN = zeros(M,4);
% count the points in each section to calculate the mean at the end.
counterpart = zeros(1,4);
countpart3 = 0;
countpart4 = 0;
countpart5 = 0;
countpart6 = 0;

for xx = 1:FieldPoints
    % TODO
    if abs(RunData(xx).XPoint) < abs(-25.50)
        SectionsLWBechara(1,1) = SectionsLWBechara(1,1) + GraphArrayLW(xx,2);
        SectionsLNBechara(1,1) = SectionsLNBechara(1,1) + GraphArrayLN(xx,2);
%         for mm = 1:M
%             SectionN(mm,1) = SectionN(mm,1) + GraphArrayN(xx,mm);
%         end
        countpart3 = countpart3 + 1;
        counterpart(1,1) = counterpart(1,1) + 1; 
        
    elseif abs(RunData(xx).XPoint) > abs(-25.50) && abs(RunData(xx).XPoint) < abs(-35.50)
        SectionsLWBechara(1,2) = SectionsLWBechara(1,2) + GraphArrayLW(xx,2);
        SectionLNBechara(1,2) = SectionsLNBechara(1,2) + GraphArrayLN(xx,2);
%         for mm = 1:M
%             SectionN(mm,2) = SectionN(mm,2) + GraphArrayN(xx,mm);
%         end
        countpart4 = countpart4 + 1;
        counterpart(1,2) = counterpart(1,2) + 1;
        
    elseif abs(RunData(xx).XPoint) > abs(-35.50) && abs(RunData(xx).XPoint) < abs(-45.50)
        SectionsLWBechara(1,3) = SectionsLWBechara(1,3) + GraphArrayLW(xx,2);
        SectionLNBechara(1,3) = SectionsLNBechara(1,3) + GraphArrayLN(xx,2);
%         for mm = 1:M
%             SectionN(mm,3) = SectionN(mm,3) + GraphArrayN(xx,mm);
%         end
        countpart5 = countpart5 + 1;
        counterpart(1,3) = counterpart(1,3) + 1;
    else
        SectionsLWBechara(1,4) = SectionsLWBechara(1,4) + GraphArrayLW(xx,2);
        SectionLNBechara(1,4) = SectionsLNBechara(1,4) + GraphArrayLN(xx,2);
%         for mm = 1:M
%             SectionN(mm,4) = SectionN(mm,4) + GraphArrayN(xx,mm);
%         end
        countpart6 = countpart6 + 1;   
        counterpart(1,4) = counterpart(1,4) + 1;
    end
end


% Averaging everything: 

FinalSectionLWBechara = zeros(1,4);
FinalSectionLNBechara = zeros(1,4);
% Nmean = zeros(4,M);

for ii = 1:4
    %TODO
    FinalSectionLWBechara(1,ii) = SectionsLWBechara(1,ii)/counterpart(1,ii);
    FinalSectionLNBechara(1,ii) = SectionLNBechara(1,ii)/counterpart(1,ii);
    
%     for mm = 1:M
%         Nmean(ii,mm)  = SectionN(mm,ii)/counterpart(1,ii);    
%     end
    
    
end


%% For averaging fiber length through the section

% First pick a section:
section =8;

probeInitial = (section-1)*layer+1;

SectionMatrixLW = zeros(layerZ,layerY); % To store fiber lengths LW
SectionMatrixSR = zeros(layerZ,layerY); % To store shear rates
SectionMatrixVis = zeros(layerZ,layerY); % To store viscosities

if section ~= sections
    SectionProbes = layer; % No of probes in a given 'x'section

    for ii = 1:SectionProbes
        %TODO
        yy = ceil(ii/layerZ);
        
        if rem(ii,layerZ) == 0
            zz = layerZ;
        else
            zz = rem(ii,layerZ);
        end
        
        SectionMatrixLW(zz,yy) = GraphArrayLW(probeInitial-1+ii,2)*1000;
        SectionMatrixSR(zz,yy) = GraphArraySR(probeInitial-1+ii,1);
        SectionMatrixVis(zz,yy) = GraphArrayVis(probeInitial-1+ii,1);
        
       
    end
else
    SectionProbes = probes - (layer*(sections-1));
end

% First invert the array: To invert matrix so that (1,1) starts at bottom left.
SectionMatrixLW = flip(SectionMatrixLW); 
SectionMatrixSR = flip(SectionMatrixSR);
SectionMatrixVis = flip(SectionMatrixVis);

% Now we plot in the section
% Plot through Z:
thicknessZLW = zeros(1,layerZ);
thicknessZSR = zeros(1,layerZ);
thicknessZVis = zeros(1,layerZ);
for zz = 1:layerZ
    thicknessZLW(1,zz) = mean(SectionMatrixLW(zz,:));
    thicknessZSR(1,zz) = mean(SectionMatrixSR(zz,:));
    thicknessZVis(1,zz) = mean(SectionMatrixVis(zz,:));
end
figure()
yyaxis left
plot(thicknessZLW)
ylabel('fiber length')
hold on
yyaxis right
plot(thicknessZSR)
ylabel('shear rate')
xlabel('thickness Z')
title('fiber length vs thickness')

% Plot through Y:
thicknessYLW = zeros(1,layerY);
thicknessYSR = zeros(1,layerY);
thicknessYVis = zeros(1,layerY);
for yy = 1:layerY
    thicknessYLW(1,yy) = mean(SectionMatrixLW(:,yy));
    thicknessYSR(1,yy) = mean(SectionMatrixSR(:,yy));
    thicknessYVis(1,yy) = mean(SectionMatrixVis(:,yy));
end
figure()
yyaxis left
plot(thicknessYLW)
ylabel('fiber length')
hold on
yyaxis right
plot(thicknessYSR)
ylabel('shear rate')
xlabel('thickness Y')
title('fiber length vs thickness')



