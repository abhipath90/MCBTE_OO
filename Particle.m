classdef Particle
    % This is the Monte Carlo particle
    %   It contains all the information a particle need to know its state,
    %   its interactions with other objects
    
    properties
        eEff     % the effective energy of the particle in simulation
        pSign    % sign of the particle
        vel      % Velocity vector
        matID    % Current material in which particle resides
        omega    % Frequency of the particle
        mode     % the mode number
        
        startPoint (3,1) double {mustBeNonNan}   % defines the starting point of current flight
        endPoint   (3,1) double {mustBeNonNan}   % defines the end point of the current flight
        t0          % Start time of the current flight
        tNext       % time to next scattering 
        
        tNext3ph    % time to next 3 phonon process
        tNextImp    % time to next impurity scattering
        isAlive     % Boolean to check if the particle is alive
        relaxCount % Number of relaxation events since birth
        
        % Don't expect more than 10000
        % scattering events. May need to find better way to initialize it
        % based on max_scatter value
        %trajAll = zeros(10000,6)      % contains starting and ending points of all the flights
        scattHist = zeros(10000,3)    % contains scattering history scatt type, if bnd bndID
        % 1 = relaxation; 2 = impurity, 3 = bnd
    end
    
    methods
        
        % For now we will use the default constructor
        %         function obj = Particle(Start_point,effectiveE,sign_of,velocity,...
        %                 start_time)
        %             % this is the constructor of the particle that is called at the
        %             % birth of the particle
        %             obj.Start = Start_point;
        %             obj.End = [];
        %             obj.Eeff = effectiveE;
        %             obj.sign = sign_of;
        %             obj.V = velocity;
        %             obj.t0 = start_time;
        %             obj.tNext = [];
        %             obj.tNext3ph = [];
        %             obj.tNextImp = [];
        %             obj.IsAlive = true;
        %         end
        
        %% Methods to give birth to the particle
        
        %%  Methods below assign values not assignes at birth
        function obj = set.tNext3ph(obj,timeToRelax)
            % Assigns the three phonon relaxation time
            obj.tNext3ph = timeToRelax;
        end
        
        function obj = set.tNextImp(obj,timeToImpure)
            % Assigns the impurity relaxation time
            obj.tNextImp = timeToImpure;
        end
        
        %% The modifier functions
        
        function obj = set_next(obj,tEnd,locEnd)
            % Assigns/modifies the flight end/start and end
            % This method is used to assign end of trajectory right after
            % the birth of particle, corresponding to the first scattering
            % event. Or it assigns the new flight to the particle.
            if(isempty(obj.tNext))
                obj.End = locEnd;
                obj.tNext = tEnd;
            else
                obj.Start = obj.End;
                obj.End = locEnd;
                obj.tNext = tEnd;
            end
        end
            
        %% Below are the functions to identify Particle's interaction with
        % Other geometric entities such as regions and boundaries
        
        function [isInside, insideID] = Check_inside(obj,Regions,time_stamp)
            % This method looks at all the Region objects and finds where
            % the particle belongs at a particular time. This method is
            % useful for transient cases where we want to know the location
            % of the particle to count its contribution to temperature and
            % flux.
            isInside = false;
            Position = obj.startPoint+ obj.vel * (time_stamp - obj.t0);
            % Checking through all the regions (currently it is only
            % defined 
            for regionID = 1:length(Regions)
                isInside = Regions(regionID).Check_inside(Regions(regionID),Position);
                if isInside
                    break;
                end
            end
            insideID = regionID;
        end
        
        
        %% Scattering function
        function obj = Relax_scatter(obj,matProp)
            % recording values for debug before scattering
            % finding where to record, only need first column to check
            % indexScat = find(~obj.scattHist(:,1) && ~obj.scattHist(:,2),1);
            indexScat = find(~obj.scattHist(:,1),1);
            obj.scattHist(indexScat,:) = [1,0,obj.matID];
            %obj.trajAll(indexScat,:) = [obj.startPoint' obj.endPoint'];
            
            % 3 phonon scattering implemented here
            obj.mode = Utils.Select_mode(matProp.cumulColl,matProp.Nmodes);
            obj.omega = matProp.freq(obj.mode);
            obj.vel = matProp.vel(obj.mode)*Utils.Draw_random('sphere');
            obj.startPoint = obj.endPoint;
            obj.t0 = obj.tNext;
            obj.tNext3ph = obj.tNext - matProp.tau(obj.mode)*log(1-rand());
            obj.tNextImp = obj.tNext - matProp.tauImp(obj.mode)*log(1-rand());
            obj.relaxCount = obj.relaxCount +1;
        end
        
        function obj = Imp_scatter(obj,matProp)
            % recording values for debug before scattering
            % finding where to record
            indexScat = find(~obj.scattHist(:,1) ,1);
            obj.scattHist(indexScat,:) = [2,0,obj.matID];
            %obj.trajAll(indexScat,:) = [obj.startPoint' obj.endPoint'];
            
            % 2 phonon/ impurity scattering is implemented
            obj.vel = vecnorm(obj.vel,2,1)*Utils.Draw_random('sphere');
            obj.startPoint = obj.endPoint;
            obj.t0 = obj.tNext;
            obj.tNextImp = obj.tNext - matProp.tauImp(obj.mode)*log(1-rand());
        end
        
        function [obj, GeoObj] = Geo_scatter(obj,GeoObj,bndID,matAllProp)
            % whenever flux is recorded across interface, Interface object
            % changes, that changes Geometry object too.
            % recording values for debug before scattering
            % finding where to record
            indexScat = find(~obj.scattHist(:,1) ,1);
            obj.scattHist(indexScat,:) = [3,bndID,obj.matID];
            %obj.trajAll(indexScat,:) = [obj.startPoint' obj.endPoint'];
            
            for ii=1:GeoObj.numBnd
                if(ii<=GeoObj.numIso)
                    % checking interaction with isothermal walls
                    if(GeoObj.Iso(ii).bndID == bndID)
                        %obj = GeoObj.Iso(ii).Scatter(obj);
                        obj.isAlive = false;
                        return;
                    end
                elseif(ii<=(GeoObj.numIso+GeoObj.numAdi))
                    % checking interaction with adiabatic walls
                    subtract = GeoObj.numIso;
                    if(GeoObj.Adi(ii-subtract).bndID == bndID)
                        obj = GeoObj.Adi(ii-subtract).Scatter(obj);
                        obj.startPoint = obj.endPoint;
                        obj.t0 = obj.tNext;
                        return;
                    end
                elseif(ii<=(GeoObj.numIso+GeoObj.numAdi+GeoObj.numPeri))
                    % checking interaction with periodic walls
                    subtract = GeoObj.numIso+GeoObj.numAdi;
                    if(GeoObj.Peri(ii-subtract).bndID == bndID)
                        obj = GeoObj.Peri(ii-subtract).Scatter(obj);
                        obj.t0 = obj.tNext;
                        return;
                    end                    
                else
                    subtract = GeoObj.numIso + GeoObj.numAdi + GeoObj.numPeri;
                    if(GeoObj.Inter(ii-subtract).bndID == bndID)
                        [obj, updatedInterface] = GeoObj.Inter(ii-subtract).Scatter(obj,matAllProp);
                        GeoObj.Inter(ii-subtract) = updatedInterface;
                        % now geometry object is updated. So need to return
                        % that.
                        obj.t0 = obj.tNext;
                        return;
                    end
                    
                end
            end
        end
    end
end
