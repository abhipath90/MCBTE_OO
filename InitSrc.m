classdef InitSrc
    % InitSrc stores sources of the type initial conditions
    %   This source will only be used in the transient case
    
    properties
        regionData          % region array where the temperature is defined
        delT                % delta temperature that is assigned
        matObj              % material properties to be used
        devELoc             % location wise deviational energy
        cumDevELoc          % absolute value cumulative energy array based on location
        devETotal           % absolute values total deviational energy of the source
        
    end
    
    methods
        function obj = InitSrc(initData,matProp)
            % InitSrc Construct an instance of this class
            %   Detailed explanation goes here
            obj.regionData= initData(:,1:6);
            obj.delT = initData(:,7);
            obj.matObj = matProp;
            
            vol = abs((initData(:,1)-initData(:,2)).*(initData(:,3)-initData(:,4)).*(initData(:,5)-initData(:,6)));
            obj.devELoc = obj.delT.*vol.*matProp.cumulBase(end);
            obj.cumDevELoc = cumsum(abs(obj.devELoc));
            obj.devETotal = obj.cumDevELoc(end);
        end
        
        function part = emit(obj)
            part = Particle();  % default constructor: all nulls
            
            % find the location/region form which the particle will be
            % originated
            rLoc = rand();
            for ii=1:length(obj.delT)
                if(rLoc<obj.cumDevELoc(ii)/obj.devETotal)
                    break;
                end
            end
            
            % choosing random point in that region/location
            x0 = obj.regionData(ii,1) + rand()*(obj.regionData(ii,2)-obj.regionData(:,1));
            y0 = obj.regionData(ii,3) + rand()*(obj.regionData(ii,4)-obj.regionData(:,3));
            z0 = obj.regionData(ii,5) + rand()*(obj.regionData(ii,6)-obj.regionData(:,5));
            
            part.pSign = sign(obj.delT(ii));
            part.matID = obj.matObj.matID;
            part.mode = Utils.Select_mode(obj.matObj.cumulBase,obj.matObj.Nmodes);
            part.omega = obj.matObj.freq(part.mode);
            
            % Finding velocity vector 
            speed = obj.matObj.vel(part.mode);
            % Finding random direction on a hemisphere above x=0 plane
            randDir = Utils.Draw_random('sphere');
            
            part.vel = speed*randDir;
             
            % Finding start point
            part.startPoint = [x0;y0;z0];
        end
    end
end

