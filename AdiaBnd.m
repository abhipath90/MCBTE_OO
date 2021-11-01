classdef AdiaBnd < Boundary
    % AdiaBnd is adiabatic boundary
    %   This class contains information about specularity of the boundary
    %   alongside the regular boundary properties
    
    properties
        Spec     % Specularity of the adiabatic boundary
    end
    
    methods
        %% Constructor
        function obj = AdiaBnd(ID,matID,Point1,Point2,Point3,Point4,S)
            if nargin < 6, S = []; end
            if nargin < 5, Point4 = []; end
            if nargin < 4, Point3 = []; end
            if nargin < 3, Point2 = []; end
            if nargin < 2, Point1 = []; end
            
            obj@Boundary(ID,matID,Point1,Point2,Point3,Point4);
            obj.Spec = S;
        end
        

        function part = Scatter(obj,part)
            rRefl = rand();
            if(rRefl<obj.Spec)
                % specular reflection
                % calling reflection routine
                part.vel = Utils.Reflection(part.vel,obj.normal);
            else
                % diffuse reflection
                % Finding random direction in the hemisphere
                randDir = Utils.Draw_random('hemisphere');
                
                % Orienting to reflect from the wall: Rotation matrix
                rotMat = Utils.Orient([1;0;0],obj.normal);
                % assigning the random velocity
                part.vel = vecnorm(part.vel,2,1)*rotMat*randDir;
            end
        end
    end
end

