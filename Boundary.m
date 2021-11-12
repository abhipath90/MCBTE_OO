classdef Boundary
    % Boundary is a geometric plane defining boundaries in the domain
    %   This will contain plane's definition alongside methods to calculate
    %   geometric properties of plane and its interaction with particle
    %   flight. It is assumed that the four points are listed
    %   counter-clockwise when viewd from the inside of the domain except
    %   when it is an interface between 2 materials
    
    properties
        bndID
        matID    % material ID towards the normal direction
        vertex1 (3,1) double {mustBeNonNan}
        vertex2 (3,1) double {mustBeNonNan}
        vertex3 (3,1) double {mustBeNonNan}
        vertex4 (3,1) double {mustBeNonNan}
        area
        normal
    end
    
    methods
        %% Constructor
        function obj = Boundary(ID,materialID,Point1,Point2,Point3,Point4)
            if nargin < 5, Point4 = []; end
            if nargin < 4, Point3 = []; end
            if nargin < 3, Point2 = []; end
            if nargin < 2, Point1 = []; end
            
            obj.bndID = ID;
            obj.matID = materialID;
            obj.vertex1 = Point1;
            obj.vertex2 = Point2;
            obj.vertex3 = Point3;
            obj.vertex4 = Point4;
            
            % Finding two sides of the rectangle as vectors
            vec1 = obj.vertex2 - obj.vertex1;
            vec2 = obj.vertex4 - obj.vertex1;
            
            % Area as cross product
            areaVec = cross(vec1,vec2,1);       % column vectors
            obj.area = vecnorm(areaVec,2,1);    % column vectors 2-norm
            obj.normal = areaVec./obj.area;
        end
            
        
        
        function randPoint = Random_point(obj)
            % this function calculates a random point on the boundary
            % Right now its just a rectangular boundary so we will be using
            % simple method of finding two random points on the two
            % diagonals and then taking an average of them
            r1 = rand();
            r2 = rand();
            randPoint = (obj.vertex1 + r1*(obj.vertex3-obj.vertex1) ... 
                + obj.vertex2 + r2*(obj.vertex4-obj.vertex2))/2;
        end
        
        
        
        function [fraction,point,ID] = interaction(obj,part)
            % calculates interaction between boundary and particle
            % Interface boundary has another function that overrides this
            % one
            % PeriBnd also overrides this
            fraction = -1; ID = obj.bndID; point=[NaN;NaN;NaN];
            indexScat = find(~part.scattHist(:,1) ,1);
            if(indexScat ~=1) % not the first scattering
                % check if the particle just scattered from this boundary
                if(part.scattHist(indexScat-1,2) == ID)
                    % just scattered from this boundary, ignore
                    return;
                end
            end
            A = part.startPoint;
            B = part.endPoint;
            
            if(dot(B-A,obj.normal)>0)
                % particle is going away from the boundary
                return;
            end
            
            if(dot(obj.normal,B-A)==0)
                % segment is parallel to the boundary
                return;
            else
                % check if the point A lies on the boundary being checked
                % for. Sometimes if the flight is too short, a repeated
                % interaction can be detected
                %********** NEED A BETTER CHECK***********************
                checkVec = obj.vertex1 - A;
                checkVec = checkVec/norm(checkVec); % for proper scaling, unit vector
                if(abs(dot(obj.normal,checkVec))< 1e-5)
                    return;
                end
                
                t = dot(obj.normal,obj.vertex1-A)/dot(obj.normal,B-A);
                
                %if(t<=eps || t>1)
                if(t<0 || t>1) % with exclusion of prev boundary 0 is possible
                    % interaction happens outside the segment
                    return;
                else
                    % I lies in the plane of the rectangle
                    I = A + t*(B-A); % interaction point
                    
                    % Lengths of the two sides of the  rectangle
                    l1 = vecnorm(obj.vertex2-obj.vertex1,2,1);
                    l2 = vecnorm(obj.vertex4-obj.vertex1,2,1);
                    
                    % check if I lies inside the rectangle
                    alpha = dot(I-obj.vertex1,obj.vertex2-obj.vertex1)/l1;
                    beta = dot(I-obj.vertex1,obj.vertex4-obj.vertex1)/l2;
                    
                    if(alpha<0 || alpha>l1 || beta<0 || beta>l2)
                        % intersection is outside rectangle
                        return;
                    else
                        fraction = t;
                        point = I;
                    end
                end
            end
        end
       


    end
end

