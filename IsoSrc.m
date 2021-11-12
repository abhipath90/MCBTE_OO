classdef IsoSrc
    % IsoSrc is the isothermal boundary
    %   It is used to emit particle as source. For each isothermal boundary
    %   defined in the boundaries, one isothermal source must be defined
    
    properties
        Temp
        Teq
        tMax
        bndObj
        matObj
        cumDevEModal        % signed modal deviational energy to pick mode
        devETotal           % Signed total deviational energy
    end
    
    methods
        function obj = IsoSrc(bnd,matProp,eqT,maxTime)
            %UNTITLED5 Construct an instance of this class
            %   Detailed explanation goes here
            obj.Teq = eqT;
            obj.bndObj = bnd;
            obj.Temp = bnd.Temp;
            obj.matObj = matProp;
            obj.tMax = maxTime;
            
            if(isempty(obj.tMax))
                % Steady state simulation
                obj.cumDevEModal = (bnd.Temp - eqT)*bnd.area.*matProp.cumulVel/4;
                obj.devETotal = obj.cumDevEModal(end);
            else
                % Transient simulation
                obj.cumDevEModal = maxTime*(bnd.Temp - eqT)*bnd.area.*matProp.cumulVel/4;
                obj.devETotal = obj.cumDevEModal(end);
            end
            
            
        end
        
        function part= emit(obj)
            % Initializes particle with relevant properties
            %   Properties like the frequency, velocity, start position are
            %   selected. Properties those do not depend on this particular
            %   source are assigned at a level above.
            
            part = Particle();  % default constructor: all nulls
            part.pSign = sign(obj.devETotal);
            part.matID = obj.matObj.matID;
            
            part.mode = Utils.Select_mode(obj.matObj.cumulVel,obj.matObj.Nmodes);
            part.omega = obj.matObj.freq(part.mode);
            
            % Finding velocity vector 
            speed = obj.matObj.vel(part.mode);
            % Finding random direction on a hemisphere above x=0 plane
            randDir = Utils.Draw_random('hemisphere');
            
            % Orienting random random vector based on boundary normal
            rotMat = Utils.Orient([1;0;0],obj.bndObj.normal);
            part.vel = speed*rotMat*randDir;
             
            % Finding start point
            part.startPoint = obj.bndObj.Random_point();
            %
        end
    end
end

