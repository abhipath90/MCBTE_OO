classdef InputData
    % Stores all the input data in raw form
    %   Later on all the transformations will be included and checks and
    %   balances implemented.
    
    properties
        matData         % A cell arrary containing material data
        bnd             % Array containing all the boundaries
        bndProp         % Array of all boundary properties in same order
        measureReg      % Array of all the detector data
        thermalGrad     % Array of thermal gradient information
        measureTime     % Vector of all the measurement times
        param           % Vector of simulations parameters
        init            % A cell array of all the initial conditions
        interface       % cell array of all the interfacial conditions
        periPair        % array containing all periodic boundary pairs
    end
    
    methods
        function obj = InputData(mat1,mat2,bndFile,propFile,regionFile,...
                gradFile,timeFile,paramFile,init1,init2,File_interface,File_peri_pair)
            % This is to store the data coming from all the files
            
            obj.matData{1,1} = load(mat1);
            obj.matData{2,1} = load(mat2);
            obj.bnd = load(bndFile);
            obj.bndProp = load(propFile);
            obj.measureReg = load(regionFile);
            obj.thermalGrad = load(gradFile);
            obj.measureTime = load(timeFile);
            obj.param = load(paramFile);
            obj.init{1,1} = load(init1);
            obj.init{2,1} = load(init2);
            obj.periPair = load(File_peri_pair);
            
            % getting interface conditions
            [dataDir, ~, ~] = fileparts(File_interface);
            fid = fopen(File_interface);
            tline = fgetl(fid); % first line discard
            numInter = str2num(fgetl(fid));  % stores number of interfaces in the simulation
            if(numInter==0)
                obj.interface = {numInter};
            else         
                tline = fgetl(fid); %discard this one too
                interBndID = str2num(fgetl(fid));  % stores boundary ids of all the interface
                % next 4 lines are comments
                for disCount=1:4
                    tline = fgetl(fid);% discarding comment lines
                end
                interDefFiles = [];
                for interCount = 1:numInter
                    interDefFiles{interCount} =  split(fgetl(fid),',');
                end
                fclose(fid);
                % Need to read all the interDefFiles here only and collect
                % data in the memory
                InterDefData = [];
                for interCount = 1:numInter
                    numFiles = length(interDefFiles{interCount});
                    for fileCount = 1:numFiles
                        InterDefData{interCount}{fileCount} = regexp(fileread([dataDir '/' strtrim(interDefFiles{interCount}{fileCount})]), '\r\n|\r|\n', 'split');
                        % strtrim removes leading whitespaces from an
                        % string
                    end
                end
                obj.interface = {numInter; interBndID; InterDefData};
            end
        end
      
        
    end
end

