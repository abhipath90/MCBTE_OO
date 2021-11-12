classdef GradSrc
    % Gradsrc is the class to store thermal gradient source
    %   This stores data and functions related to thermal gradient source.
    %   This includes the deviational energy and method to produce
    %   particles
    
    properties
        grad        % a (3,1) vector to store gradient vector
        %vol         % volume of the domain it is applied in
        matObj      % material object for the material of the domain
        regionObj   % region object to specify the region
        tMax
        cumDevEModal
        devETotal
    end
    
    methods
        function obj = GradSrc(matProp,gradient,region,maxT)
            %GradSrc Construct an instance of this class
            %   Detailed explanation goes here
            
            obj.grad= gradient;     % gradient has to be (3,1) vector
            obj.regionObj = Region(region(1),region(2),region(3),region(4),region(5),region(6),matProp.matID);
    %        obj.vol = volume;
            obj.matObj = matProp;
            obj.tMax = maxT;
            volume = obj.regionObj.volume;
            
            if(isempty(maxT))
                % Steady state simulation
                obj.cumDevEModal = vecnorm(gradient,2,1)*volume.*matProp.cumulVel/2;
                obj.devETotal = obj.cumDevEModal(end);
            else
                % transient simulation
                obj.cumDevEModal = maxT*vecnorm(gradient,2,1)*volume.*matProp.cumulVel/2;
                obj.devETotal = obj.cumDevEModal(end);
            end
        end
        
        function part = emit(obj)
            
            part = Particle();  % default constructor: all nulls
            part.pSign = sign(rand()-0.5); % equal probability of each sign
            part.matID = obj.matObj.matID;

            part.mode = Utils.Select_mode(obj.matObj.cumulVel,obj.matObj.Nmodes);
            part.omega = obj.matObj.freq(part.mode);
                        
            % Finding velocity vector 
            speed = obj.matObj.vel(part.mode);
            % Finding random direction on a hemisphere above x=0 plane
            randDir = Utils.Draw_random('hemisphere');
            
            % Orienting random random vector based on unit vector in
            % gradient direction
            rotMat = Utils.Orient([1;0;0],obj.grad/vecnorm(obj.grad,2,1));
            
            % particle carries energy againnt the thermal gradient that
            % effect is achieved by multiplying with -psign
            part.vel = -part.pSign*speed*rotMat*randDir;
             
            % Finding start point
            part.startPoint = obj.regionObj.Random_point();
            %
        end
    end
end

