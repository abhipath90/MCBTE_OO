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
    end
    
    methods
        function obj = InputData(mat1,mat2,bndFile,propFile,regionFile,...
                gradFile,timeFile,paramFile,init1,init2,File_interface)
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
            
            % getting interface conditions
            fid = fopen(File_interface);
            tline = fgetl(fid); % first line discard
            numInter = str2num(fgetl(fid));
            if(numInter==0)
                obj.interface = {numInter};
            else         
                tline = fgetl(fid); %discard this one too
                interBndID = str2num(fgetl(fid));
                tline = fgetl(fid);% discarding again
                interDataNum = str2num(fgetl(fid));
                % get new fid
                fclose(fid);
                %fid = fopen(File_interface);
                %interData = cell2mat(textscan(fid,'%f %f %f %f','HeaderLines',6));
                interData = readmatrix(File_interface,'NumHeaderLines',6,'Delimiter',',');
                obj.interface = {numInter; interBndID;interDataNum;interData};
            end
        end
        
        
    end
end

