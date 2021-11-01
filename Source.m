classdef Source
    % This class contains the source definition
    %   One class to keep track of all the sources and it will choose which
    %   source to be chosen to emit particle
    
    properties
        numSrc
        numIso
        numInit
        numGrad
        isoSrc      % Array of isothermal sources
        initSrc     % Array of initial conditions
        gradSrc     % Array of thermal gradient sources
        tMax
        simType
        cumDevE      % cumulative deviational energy for all the sources
        devE        % total deviational energy in the simulation
    end
    
    methods
        function obj = Source(RawData,simObj,GeoObj,MaterialInfo)
            % Create all sources and store related information
            %   This method will read all the files related to sources and
            %   then create all the source.
            eqT = simObj.Teq;
            if(strcmpi(simObj.typeSim,'Steady'))
                obj.tMax = [];
                obj.simType = 'Steady';
            else
                obj.tMax = simObj.tMax;
                obj.simType = 'Transient';
            end
            
            % creating isothermal sources
            obj.numIso = GeoObj.numIso;
            obj.isoSrc = [];
            if(obj.numIso>0)
                for ii=1:obj.numIso
                    bndObj = GeoObj.Iso(ii);
                    matID = bndObj.matID;
                    mat = MaterialInfo.matProp(matID);
                    obj.isoSrc = [obj.isoSrc; IsoSrc(bndObj,mat,eqT,obj.tMax)];
                end
            end
            
            % creating thermal gradient sources
            % currently only supports one material and one domain whoose
            % volume must be given in the SimParam.txt file.
            obj.numGrad = 0;
            obj.gradSrc = [];
            if(~isempty(RawData.thermalGrad))
                [obj.numGrad,~] = size(RawData.thermalGrad);
                for ii=1:obj.numGrad
                    % Taking default material, the first one. May need to
                    % be generalized later on.
                    grad = (RawData.thermalGrad(ii,3:5))';      % gradient value
                    gradRegion = RawData.thermalGrad(ii,6:11);  % region
                    gradMat = RawData.thermalGrad(ii,12);       % material ID
                    obj.gradSrc = [obj.gradSrc; GradSrc(MaterialInfo.matProp(gradMat),grad,gradRegion,obj.tMax)];
                end
            end
            
            % creating intial conditions sources
            
            obj.numInit = 0;
            obj.initSrc = [];
            
            % right now there is one initial condition file for each
            % material arranged in the order of their material ID. Later a
            % more robust way can be added
            
            % finding number of defined initial coditions and their matIDs
            matInitIDs = find(~cellfun(@isempty,RawData.init));
            
            if(~isempty(matInitIDs))
                obj.numInit = length(matInitIDs);
                for ii=1:length(matInitIDs)
                    id = matInitIDs(ii);
                    initDataIn = cell2mat(RawData.init(id));
                    obj.initSrc = [obj.initSrc; InitSrc(initDataIn,MaterialInfo.matProp(id))];
                end
            end
            
            % All sources
            obj.numSrc = obj.numIso + obj.numInit + obj.numGrad;
            
            % cumulative deviational energy for choosing source to emit the
            % order will be iso then grad and then init
            obj.cumDevE = zeros(obj.numSrc,1);
            count=1;
            % putting isothermal sources
            for ii=1:obj.numIso
                if count==1
                    obj.cumDevE(count) = abs(obj.isoSrc(ii).devETotal);
                else
                    obj.cumDevE(count) = obj.cumDevE(count-1) + abs(obj.isoSrc(ii).devETotal);
                end
                count = count +1;
            end
            
            % putting gradient source
            for ii=1:obj.numGrad
                if count == 1
                    obj.cumDevE(count) = abs(obj.gradSrc(ii).devETotal);
                else
                    obj.cumDevE(count) = obj.cumDevE(count-1) + abs(obj.gradSrc(ii).devETotal);
                end
                count = count +1;
            end
            
            % putting intial conditions source
            for ii=1:obj.numInit
                if count ==1
                    obj.cumDevE(count) = abs(obj.initSrc(ii).devETotal);
                else
                    obj.cumDevE(count) = obj.cumDevE(count-1) + abs(obj.initSrc(ii).devETotal);
                end
                count = count +1;
            end
            
            obj.devE = obj.cumDevE(end);
        end
        
        %% Emission of particle from source
        function particle = Emit(obj,N)
            % This method will create a particle based on the input
            % properties
            %   Method chooses the correct source to emit particle from and
            %   then emits the particle from that by assigning all
            %   appropriate initial conditions
            
            % Choose the source to emit particle
            cumProb = abs(obj.cumDevE)./obj.devE;
            rSrc = rand();
            for srcIndex=1:length(cumProb)
                if rSrc < cumProb(srcIndex)
                    break;
                end
            end
            % Assumption: the sources are arranged in the order of
            % isothermal, gradient and initial conditions
            if(srcIndex<=obj.numIso)
                % emit from isothermal boundaries
                isoIndex = srcIndex;
                particle = obj.isoSrc(isoIndex).emit();
                
                % Assigning properties not assgined by the isothermal
                % source
                particle.eEff = obj.devE/N;
                if(strcmpi(obj.simType,'Steady'))
                    particle.t0 = 0;
                else
                    particle.t0 = obj.tMax*rand();
                end
                    
                particle.tNext3ph = particle.t0 -obj.isoSrc(isoIndex).matObj.tau(particle.mode)*log(1-rand());
                particle.tNextImp = particle.t0 -obj.isoSrc(isoIndex).matObj.tauImp(particle.mode)*log(1-rand());
                
                particle.isAlive = true;
                particle.relaxCount = 0;
                
            elseif(srcIndex<=(obj.numIso+obj.numGrad))
                % emit from gradient source
                gradIndex = srcIndex - obj.numIso;
                particle = obj.gradSrc(gradIndex).emit();
                
                % Assigning properties not assgined by the gradient
                % source
                particle.eEff = obj.devE/N;
                if(strcmpi(obj.simType,'Steady'))
                    particle.t0 = 0;
                else
                    particle.t0 = obj.tMax*rand();
                end
                
                particle.tNext3ph = particle.t0 -obj.gradSrc(gradIndex).matObj.tau(particle.mode)*log(1-rand());
                particle.tNextImp = particle.t0 -obj.gradSrc(gradIndex).matObj.tauImp(particle.mode)*log(1-rand());
                particle.isAlive = true;
                particle.relaxCount = 0;
                
            else
                % emit from initial conditions
                initIndex = srcIndex - (obj.numIso + obj.numGrad);
                particle = obj.initSrc(initIndex).emit();
                
                % Assigning properties not assgined by the gradient
                % source
                particle.eEff = obj.devE/N;
                if(strcmpi(obj.simType,'Steady'))
                    error('There should not be an intial condition for steady state simulation\n');
                else
                    % If it originates from initial conditions, it will
                    % always have t0 = 0
                    particle.t0 = 0;
                end
                
                particle.tNext3ph = particle.t0 -obj.initSrc(initIndex).matObj.tau(particle.mode)*log(1-rand());
                particle.tNextImp = particle.t0 -obj.initSrc(initIndex).matObj.tauImp(particle.mode)*log(1-rand());
                particle.isAlive = true;
                particle.relaxCount = 0;
            end
            
        end
    end
end

