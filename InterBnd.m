classdef InterBnd < Boundary
    % This is interface boundary between two materials
    %   It is an internal boundary in the geometric domain and will be used
    %   to store interfaces between simmilar or dissimmilar materials. It
    %   will also store methods to interact with particles such as
    %   calculating if interaction happens during the flight and what
    %   happens to the particle incident on the interface
    
    properties
        %matID1          % material towards the normal direction
        matID2          % material on the other side
        T12             % transmission coefficient from mat1 to mat2
        T21             % transmission coefficient from mat2 to mat1
        S12             % specularity of transmission and reflection from mat1
        S21             % specularity of transmission and reflection from mat2
        %intType         % interface type. To determine characteristics
    end
    
    methods
        function obj = InterBnd(ID,matID,Point1,Point2,Point3,Point4,mat2,type,matObj)
            % Constructor for interface boundary object
            obj@Boundary(ID,matID,Point1,Point2,Point3,Point4);
            
            %obj.matID1 = mat1;
            % This complicated line is here just to use matObj once.
            obj.matID2 = matObj.matProp(mat2).matID;%mat2;
            %obj.intType = type;
            
            % These variables may have some trailing NaNs if the number of
            % data points in mat1 and mat2 are not same. However they
            % should never be accessed anywhere in the problem if
            % everything works as it is supposed to work.
            % rmmissing removes any NaNs from the vectors.
            obj.T12 = rmmissing(type(:,1));
            obj.T21 = rmmissing(type(:,2));
            obj.S12 = rmmissing(type(:,3));
            obj.S21 = rmmissing(type(:,4));
        end
        
        
        
        function [fraction,point,ID] = interaction(obj,part)
            % this function is redfined here becuase interface will only be
            % defined once but particle can approach it from both the sides
            fraction = -1; ID = obj.bndID; point=[NaN;NaN;NaN];
                        
            A = part.startPoint;
            B = part.endPoint;
            
            % commenting the approach check part here.
            
%             if(dot(B-A,obj.normal)>0)
%                 % particle is going away from the boundary
%                 return;
%             end
            
            if(dot(obj.normal,B-A)==0)
                % segment is parallel to the boundary
                return;
            else
                t = dot(obj.normal,obj.vertex1-A)/dot(obj.normal,B-A);
                
                if(t<=eps || t>1)
                    % interaction happens outside the segment
                    return;
                else
                    I = A + t*(B-A); % interaction point
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
        
        
        function part = Scatter(obj,part,matAllProp)
            if(part.matID == obj.matID)
                %normalToUse = obj.normal;
                % particle approaching from mat1 side
                if(size(obj.T12)>1)
                    % modal values are available
                    transProb = obj.T12(part.mode);
                else
                    transProb = obj.T12;
                end
                if(size(obj.S12)>1)
                    % modal values are available
                    spec = obj.S12(part.mode);
                else
                    spec = obj.S12;
                end
            elseif(part.matID == obj.matID2)
                %normalToUse = -obj.normal;
                % particle approaching from mat2 side
                if(size(obj.T21)>1)
                    % modal values are available
                    transProb = obj.T21(part.mode);
                else
                    transProb = obj.T21;
                end
                if(size(obj.S21)>1)
                    % modal values are available
                    spec = obj.S21(part.mode);
                else
                    spec = obj.S21;
                end
            else
                error('Particle and interface mat IDs do not match');                
            end
            
            rTrans = rand();
            rSpec = rand();
            
            if(rTrans<transProb)
                if(rSpec<spec)
                    % specular transimission
                    partDir = part.vel/vecnorm(part.vel,2,1);
                    
%                     if(obj.matID == obj.matID2)
%                         % both side materials are same. No need to resample
%                         % the properties
%                         part.startPoint = part.endPoint;    
                    if(part.matID == obj.matID)
                        % particle going into material 2
                        part.matID = obj.matID2;
                        part.mode = dsearchn(matAllProp(part.matID).freq,part.omega);
                        part.omega = matAllProp(part.matID).freq(part.mode);
                        part.vel = matAllProp(part.matID).vel(part.mode)*partDir;
                        part.startPoint = part.endPoint;
                        % particle will be drawn new in new material
                        part.tNext3ph = part.tNext - matAllProp(part.matID).tau(part.mode)*log(1-rand());
                        part.tNextImp = part.tNext - matAllProp(part.matID).tauImp(part.mode)*log(1-rand());
                        part.relaxCount = part.relaxCount +1;
                    else
                        % particle going into material 1
                        part.matID = obj.matID;
                        part.mode = dsearchn(matAllProp(part.matID).freq,part.omega);
                        part.omega = matAllProp(part.matID).freq(part.mode);
                        part.vel = matAllProp(part.matID).vel(part.mode)*partDir;
                        part.startPoint = part.endPoint;
                        % particle will be drawn new in new material
                        part.tNext3ph = part.tNext - matAllProp(part.matID).tau(part.mode)*log(1-rand());
                        part.tNextImp = part.tNext - matAllProp(part.matID).tauImp(part.mode)*log(1-rand());
                        part.relaxCount = part.relaxCount +1;
                        
                    end
                else
                    % diffuse transmission
                    % new particle direction will be drawn as if it is
                    % emitted from the other side of the interface
                    % diffusively
                    if(dot(part.vel,obj.normal)>0)
                        normalToUse = obj.normal;
                    else
                        normalToUse = - obj.normal;
                    end
                    randDir = Utils.Draw_random('hemisphere');
                    rotMat = Utils.Orient([1;0;0],normalToUse);
                    partDir = rotMat*randDir;
                    
%                     if(obj.matID == obj.matID2)
%                         % The material on both sides are same. There is no
%                         % need to sample the properties of the particle
%                         part.startPoint = part.endPoint;
%                         part.vel = vecnorm(part.vel,2,1)*partDir;
                    
                    if(part.matID == obj.matID)
                        % particle going into material 2
                        part.matID = obj.matID2;
                        part.mode = dsearchn(matAllProp(part.matID).freq,part.omega);
                        part.omega = matAllProp(part.matID).freq(part.mode);
                        part.vel = matAllProp(part.matID).vel(part.mode)*partDir;
                        part.startPoint = part.endPoint;
                        % particle will be drawn new in new material
                        part.tNext3ph = part.tNext - matAllProp(part.matID).tau(part.mode)*log(1-rand());
                        part.tNextImp = part.tNext - matAllProp(part.matID).tauImp(part.mode)*log(1-rand());
                        part.relaxCount = part.relaxCount +1;
                    else
                        % particle going into material 1
                        part.matID = obj.matID;
                        part.mode = dsearchn(matAllProp(part.matID).freq,part.omega);
                        part.omega = matAllProp(part.matID).freq(part.mode);
                        part.vel = matAllProp(part.matID).vel(part.mode)*partDir;
                        part.startPoint = part.endPoint;
                        % particle will be drawn new in new material
                        part.tNext3ph = part.tNext - matAllProp(part.matID).tau(part.mode)*log(1-rand());
                        part.tNextImp = part.tNext - matAllProp(part.matID).tauImp(part.mode)*log(1-rand());
                        part.relaxCount = part.relaxCount +1;
                        
                    end
                end
            else
                % reflection of the praticle
                if(dot(part.vel,obj.normal)>0)
                    normalToUse = -obj.normal;
                else
                    normalToUse = obj.normal;
                end
                
                if(rSpec<spec)
                    % specular reflection
                    part.vel = Utils.Reflection(part.vel,normalToUse);
                    part.startPoint = part.endPoint;
                else
                    % reflected using normalToUse
                    randDir = Utils.Draw_random('hemisphere');
                    rotMat = Utils.Orient([1;0;0],normalToUse);
                    partDir = rotMat*randDir;
                    
                    part.vel = vecnorm(part.vel,2,1)*partDir;
                    part.startPoint = part.endPoint;
                end
            end
                        
        end
    end
end

