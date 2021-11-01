classdef Point
    % A generic point in 3D
    %   Detailed explanation goes here
    
    properties
        vec (3,1) {mustBeNonNan}
    end
    
    methods
        function obj = Point(a,b,c)
            obj.vec=[a;b;c];
        end
        
%         function outputArg = method1(obj,inputArg)
%             %METHOD1 Summary of this method goes here
%             %   Detailed explanation goes here
%             outputArg = obj.Property1 + inputArg;
%         end
    end
end

