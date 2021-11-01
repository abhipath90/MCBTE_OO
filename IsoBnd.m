classdef IsoBnd < Boundary
    % IsoBnd is isothermal boundary
    %   This class contains information about boundary temperature
    %   alongside the regular boundary properties
    
    properties
        Temp
    end
    
    methods
        %% Constructor
        function obj = IsoBnd(ID,matID,Point1,Point2,Point3,Point4,T)
            if nargin < 6, T = []; end
            if nargin < 5, Point4 = []; end
            if nargin < 4, Point3 = []; end
            if nargin < 3, Point2 = []; end
            if nargin < 2, Point1 = []; end
            
            obj@Boundary(ID,matID,Point1,Point2,Point3,Point4);
            obj.Temp = T;
        end


    end
end

