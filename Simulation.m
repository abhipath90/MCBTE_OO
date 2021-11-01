function [time] = Simulation(TestCaseFolder)


% Input folder
%TestCaseFolder = 'InputFiles';
outputIdentifier = ['Output' num2str(labindex)];

% setting up output and progress files
outputMessages = fopen([TestCaseFolder '/Output' num2str(labindex) '.out'], 'w');
fprintf(outputMessages,['Hello! I am ' num2str(labindex) ' of ' num2str(numlabs) ' workers \n']);
progress = fopen([TestCaseFolder '/Progress' num2str(labindex) '.txt'], 'w');


% Main simulations how it is structured
% simulation objects are created at the main worker and those are sent to
% other workers for simulation progression
%% Step 1 Setup input data object
% This object only contains raw data supplied as input to the simulation.
% Here we will add data checks later

% File names to be imported 
File_mat1 = [ TestCaseFolder '/mat_data.txt']; % contains material data
File_mat2 = [ TestCaseFolder '/mat_data2.txt'];
File_bnd = [ TestCaseFolder '/Boundary.txt']; % Only file for all the boundaries
File_bnd_prop = [ TestCaseFolder '/Boundary_prop.txt']; % Properties of all the boundaries
File_measure_reg = [ TestCaseFolder '/Measure_regions.txt']; % Regions where output is calculated
File_thermal_grad = [ TestCaseFolder '/Thermal_gradient.txt']; % Thermal gradient source
File_measure_time = [ TestCaseFolder '/Measure_times.txt']; % Measurement times for transient simulation
File_param = [ TestCaseFolder '/Sim_param.txt']; % Other simulation paramters
File_init1 = [ TestCaseFolder '/Initial_conditions1.txt'];   % in mat1
File_init2 = [ TestCaseFolder '/Initial_conditions2.txt'];   % in mat2
File_interface = [TestCaseFolder '/Interface_data.txt'];  % interface properties

if labindex == 1
    rawDataObject = InputData(File_mat1,File_mat2,File_bnd,File_bnd_prop,...
        File_measure_reg,File_thermal_grad,File_measure_time,File_param,...
        File_init1,File_init2,File_interface);
    
    
    %% Step 2 Setup Sim object
    % This object sets up number of particles, max scattering, impurtity
    % on/off, Steady/transient, Teq
    simParamObject = Sim(rawDataObject);
    
    
    %% Step 3 create materials object
    % This object store data about the material, performs necessary checks for
    % material data and calculates various distributions required by other
    % parts of the simulations
    materialInfoObj = MatInfo(rawDataObject);
    
    %% Step 4 create geometry object
    % This object stores the detailed information about all types of geometry
    % features. Boundary and Boundary types. All the rules with which particle
    % interacts with geometry are stored in this object
    geometryObject = Geometry(rawDataObject,materialInfoObj);
    
    
    %% Step 5 create Source object
    % This object keeps track of all the sources and the related information.
    % Deviational energies, cumulative deviational energy etc. It will also
    % calculate total deviational energy and effective deviational energy for
    % the particle. this object will release particle into domain.
    sourceObject = Source(rawDataObject,simParamObject,geometryObject,materialInfoObj);
    %% Step 6 setup measurement object
    % This object keeps track of all the measurements that are needed. This
    % will have all the detectors and methods that will be used to incoroporate
    % changes to their data from particle
    DetectObject = Detector(simParamObject,rawDataObject.measureReg,rawDataObject.measureTime,materialInfoObj.matProp);
    
end
%% Step 7 Send data to all the workers

if labindex == 1
    objectOfObjects = {rawDataObject;simParamObject;materialInfoObj;...
                    geometryObject;sourceObject;DetectObject};
    [~] = labBroadcast(1,objectOfObjects);
else
    objectOfObjects = labBroadcast(1);
end

% extracting data for all the objects
rawDataObject = objectOfObjects{1};
simParamObject = objectOfObjects{2};
materialInfoObj = objectOfObjects{3};
geometryObject = objectOfObjects{4};
sourceObject = objectOfObjects{5};
DetectObject = objectOfObjects{6};


%% Step 8  Particle object
% This object contains the information about the particle and how its
% dynamics is chagned through effects other than interaction with
% geometry. This includes 3 phonon processes, 2 phonon processes, keeping
% track of the material in which it is residing to compute its velocity
% etc.

%{
Step one: emit particle from the source object
Step two: while it stays alive
            move it till next scattering event
            record contribution along the flight
            update properties and count relaxation events
Step three: kill particle if it statisfies death criteria
%}
for ii=labindex:numlabs:simParamObject.N
    %disp(ii);
    if(mod(ii*1000,simParamObject.N)==0)
        fprintf(progress,'Current %f %% \n', ii*100/simParamObject.N);
        disp(num2str(ii*100/simParamObject.N));
    end
    
    part = sourceObject.Emit(simParamObject.N);
    
    while(part.isAlive)
        
        part.tNext = min(part.tNext3ph,part.tNextImp);
                
        part.endPoint = part.startPoint + part.vel*(part.tNext - part.t0);
        
        % check for interaction with geometry
        % Find minimum positive time at which it interacts with a boundary
        % and which boundary it is
        
        isBndScatter = false;
        [time,point,bndID]= geometryObject.Geometric_interaction(part);
        
        % Choose which scattering will happen first and reset end point
        if(~isempty(time) && time < part.tNext)
            % && only evaluates second argument when first argument is true
            % this solves issue with empty time
            part.tNext = time;
            %part.endPoint = part.startPoint + part.vel*(part.tNext-part.t0);
            part.endPoint = point;
            isBndScatter = true;
        end
        
        % Calculate contributions
        DetectObject = DetectObject.Record_contributions(part,materialInfoObj.matProp);
        
        % Change particle properties based on the scattering event
        
        % Choose scattering
        
        if(isBndScatter)
            % boundary scattering
            part = part.Geo_scatter(geometryObject,bndID,materialInfoObj.matProp);
        else
            if(abs(part.tNext-part.tNext3ph)<2*eps)
                % 3 phonon scattering
                part = part.Relax_scatter(materialInfoObj.matProp(part.matID));
            elseif(abs(part.tNext-part.tNextImp)<2*eps)
                % impurity scattering
                part = part.Imp_scatter(materialInfoObj.matProp(part.matID));
            else
                error('no scattering was chosen! something is not right');
            end
        end
        
        % check for maxscatt and max time
        if(strcmpi(simParamObject.typeSim,'Steady'))
            if(part.relaxCount == simParamObject.maxScatter)
                part.isAlive = false;
            end
        else
            if(part.t0 > simParamObject.tMax)
                part.isAlive = false;
            end
        end
    end
    

end



%% Step 9 Setup data to be written to the disk
DetectObject.Write_output(TestCaseFolder,outputIdentifier);


fclose(progress);
fclose(outputMessages);
end