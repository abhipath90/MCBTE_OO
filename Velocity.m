classdef Velocity
    % A generic point in 3D
    %   Detailed explanation goes here
    
    properties
        Vx
        Vy
        Vz
    end
    
    methods
        function obj = Velocity(a,b,c)
            obj.Vx=a;
            obj.Vy=b;
            obj.Vz=c;
        end
        
%         function outputArg = method1(obj,inputArg)
%             %METHOD1 Summary of this method goes here
%             %   Detailed explanation goes here
%             outputArg = obj.Property1 + inputArg;
%         end
    end
end
