classdef PeriBnd < Boundary
    % PeriBnd is periodic boundary
    %   This class contains information about translational vector
    %   alongside the regular boundary properties
    
    properties
        translate
        pairBnd
    end
    
    methods
        %% Constructor
        function obj = PeriBnd(ID,matID,Point1,Point2,Point3,Point4,translation,pairID)
            if nargin < 7, translation = []; end
            if nargin < 6, Point4 = []; end
            if nargin < 5, Point3 = []; end
            if nargin < 4, Point2 = []; end
            if nargin < 3, Point1 = []; end
            
            obj@Boundary(ID,matID,Point1,Point2,Point3,Point4);
            obj.translate = translation;
            obj.pairBnd = pairID;
        end
        
        function part = Scatter(obj,part)
            part.startPoint = part.endPoint + obj.translate;
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
                
                if(part.scattHist(indexScat-1,2) == obj.pairBnd)
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

