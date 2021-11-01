classdef PeriBnd < Boundary
    % PeriBnd is periodic boundary
    %   This class contains information about translational vector
    %   alongside the regular boundary properties
    
    properties
        translate
    end
    
    methods
        %% Constructor
        function obj = PeriBnd(ID,matID,Point1,Point2,Point3,Point4,translation)
            if nargin < 6, translation = []; end
            if nargin < 5, Point4 = []; end
            if nargin < 4, Point3 = []; end
            if nargin < 3, Point2 = []; end
            if nargin < 2, Point1 = []; end
            
            obj@Boundary(ID,matID,Point1,Point2,Point3,Point4);
            obj.translate = translation;
        end
        
        function part = Scatter(obj,part)
            part.startPoint = part.endPoint + obj.translate;
        end

    end
end

