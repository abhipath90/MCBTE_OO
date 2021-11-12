classdef InterBnd < Boundary
    % This is interface boundary between two materials
    %   It is an internal boundary in the geometric domain and will be used
    %   to store interfaces between simmilar or dissimmilar materials. It
    %   will also store methods to interact with particles such as
    %   calculating if interaction happens during the flight and what
    %   happens to the particle incident on the interface
    
    properties
        %matID1          % material towards the normal direction
        %(deprecated as this information is inherited from Boundary object)
        matID2          % material on the other side
        T12             % transmission coefficient from mat1 to mat2
        T21             % transmission coefficient from mat2 to mat1
        S12             % specularity of transmission and reflection from mat1
        S21             % specularity of transmission and reflection from mat2
        intType         % interface type. 1 for elastic 2 for inelastic
        fluxAcross      % Net flux crossing the interface in the simulation 
        % Properties below this line are only needed for inelastic
        % transmission and reflection for others these will be empty
        % variables
        TransProbs12    % Transmission probabilities from mat1 to mat2. Cell array.
        TransIndices12  % Corresponding indices of mat2. Cell array
        RefProbs12      % Reflection probabilities when going from mat1 to mat2. Cell array
        RefIndices12    % Corresponding indices of mat1. Cell array
        
        TransProbs21    % Transmission probabilities from mat2 to mat1. Cell array.
        TransIndices21  % Corresponding indices of mat1. Cell array
        RefProbs21      % Reflection probabilities when going from mat2 to mat1. Cell array
        RefIndices21    % Corresponding indices of mat2. Cell array
    end
    
    methods
        function obj = InterBnd(ID,matID,Point1,Point2,Point3,Point4,mat2,type,matObj)
            % Constructor for interface boundary object
            obj@Boundary(ID,matID,Point1,Point2,Point3,Point4);
            
            %obj.matID1 = mat1;
            % This complicated line is here just to use matObj once.
            obj.matID2 = matObj.matProp(mat2).matID;%mat2;
            
            % type contains all the data required to assign characteristics
            % to the interface
            if(length(type)>1)
                obj.intType = 2; % inelastic boundary
            else
                obj.intType = 1; % elastic boundary
            end

            obj.fluxAcross = 0; % initializing with a zero flux
            
            % Common for both elastic and inelastic boundaries
            % combination of fileread and regexp add one empty element at
            % the end of the cell array thus (length -1)
            transSpec = zeros(length(type{1})-1,4);
            for ii=1:(length(type{1})-1)
                transSpec(ii,:) = str2num(type{1}{ii});
            end
            % These variables may have some trailing NaNs if the number of
            % data points in mat1 and mat2 are not same. However they
            % should never be accessed anywhere in the problem if
            % everything works as it is supposed to work.
            % rmmissing removes any NaNs from the vectors.
            obj.T12 = rmmissing(transSpec(:,1));
            obj.T21 = rmmissing(transSpec(:,2));
            obj.S12 = rmmissing(transSpec(:,3));
            obj.S21 = rmmissing(transSpec(:,4));
            
            % if inelastic interface
            if(obj.intType == 2)
                % combination of fileread and regexp add one empty element at
                % the end of the cell array thus (end -1)
                obj.TransProbs12 = type{2}(1:end-1)';
                obj.TransIndices12 = type{3}(1:end-1)';
                obj.RefProbs12 =type{4}(1:end-1)';
                obj.RefIndices12 = type{5}(1:end-1)';
                obj.TransProbs21 = type{6}(1:end-1)';
                obj.TransIndices21 = type{7}(1:end-1)';
                obj.RefProbs21 = type{8}(1:end-1)';
                obj.RefIndices21 = type{9}(1:end-1)';
            end
        end
        
        
        
        function [fraction,point,ID] = interaction(obj,part)
            % this function is redfined here becuase interface will only be
            % defined once but particle can approach it from both the sides
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
            
            % commenting the approach check part here. As interface accepts
            % interactions from both sides. Defined normal is not useful
            
%             if(dot(B-A,obj.normal)>0)
%                 % particle is going away from the boundary
%                 return;
%             end
            
            if(dot(obj.normal,B-A)==0)
                % segment is parallel to the boundary
                return;
            else
                
%                 % check if the point A lies on the boundary being checked
%                 % for. Sometimes if the flight is too short, a repeated
%                 % interaction can be detected
%                 checkVec = obj.vertex1 - A;
%                 checkVec = checkVec/norm(checkVec); % for proper scaling, unit vector
%                 if(abs(dot(obj.normal,checkVec))< 1e-5)
%                     return;
%                 end
                
                t = dot(obj.normal,obj.vertex1-A)/dot(obj.normal,B-A);
                
                if(t<=eps || t>1)
                %if(t<0 || t>1) % 0 is allowed if prev boundary interaction is checked above
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
        
        
        function [part, obj] = Scatter(obj,part,matAllProp)
            %flux = 0;   % flux exchanged in the interface scattering.
            if(part.matID == obj.matID)
                %normalToUse = obj.normal;
                % particle approaching from mat1 side
                if(length(obj.T12)>1)
                    % modal values are available
                    transProb = obj.T12(part.mode);
                else
                    transProb = obj.T12;
                end
                if(length(obj.S12)>1)
                    % modal values are available
                    spec = obj.S12(part.mode);
                else
                    spec = obj.S12;
                end
            elseif(part.matID == obj.matID2)
                %normalToUse = -obj.normal;
                % particle approaching from mat2 side
                if(length(obj.T21)>1)
                    % modal values are available
                    transProb = obj.T21(part.mode);
                else
                    transProb = obj.T21;
                end
                if(length(obj.S21)>1)
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
                % particle transmitted flux will be recored in the end
                % material ID will be changed
                if(part.matID == obj.matID)
                    % Particle is in the side where normal points going on
                    % the other side
                    part.matID = obj.matID2;
                    incMat = 1; % useful for inelastic transport
                else
                    % Particle is in the side where normal does not point,
                    % going to the side where normal points
                    part.matID = obj.matID;
                    incMat = 2; % useful for inelastic transport
                end
                
                % Picking material properties if they change
                if ( obj.intType == 1)
                    part.mode = dsearchn(matAllProp(part.matID).freq,part.omega);
                else
                    % Find_mode(mode, incMat, isTrans)
                    part.mode = obj.Find_mode(part.mode, incMat, true );
                end
                part.omega = matAllProp(part.matID).freq(part.mode);
                speed = matAllProp(part.matID).vel(part.mode);
                % particle will be drawn new in new material
                part.tNext3ph = part.tNext - matAllProp(part.matID).tau(part.mode)*log(1-rand());
                part.tNextImp = part.tNext - matAllProp(part.matID).tauImp(part.mode)*log(1-rand());
                part.relaxCount = part.relaxCount +1;
                
                % resetting starting point
                part.startPoint = part.endPoint;
                
                % finding new direction if it changes
                % note that part vel has not been touched yet. Make
                % sure that it is indeed the case.
                if(rSpec<spec)
                    partDir = part.vel/vecnorm(part.vel,2,1);
                else
                    if(dot(part.vel,obj.normal)>0)
                        normalToUse = obj.normal;
                    else
                        normalToUse = - obj.normal;
                    end
                    randDir = Utils.Draw_random('hemisphere');
                    rotMat = Utils.Orient([1;0;0],normalToUse);
                    partDir = rotMat*randDir;
                end
                part.vel = speed*partDir;

                
                % ******************************************************************************************************************************
                % Flux transmitted
                flux = part.pSign*part.eEff*part.vel;
                
                obj.fluxAcross = obj.fluxAcross + flux;
                % ******************************************************************************************************************************
            else
                % reflection of the praticle, No flux generated across
                % material ID will remain the same
                if(part.matID == obj.matID)
                    % Particle is in the side where normal points
                    incMat = 1; % useful for inelastic transport
                else
                    % Particle is in the side where normal does not point
                    incMat = 2; % useful for inelastic transport
                end
                
                % Picking material properties if they change
                if ( obj.intType == 1)
                    % elastic reflection only the new direction needs to be
                    % chosen
                    speed = matAllProp(part.matID).vel(part.mode);
                else
                    % inelastic refelction new mode needs to be chosen
                    % Find_mode(mode, incMat, isTrans)
                    part.mode = obj.Find_mode(part.mode, incMat, false );
                    part.omega = matAllProp(part.matID).freq(part.mode);
                    speed = matAllProp(part.matID).vel(part.mode);
                    part.tNext3ph = part.tNext - matAllProp(part.matID).tau(part.mode)*log(1-rand());
                    part.tNextImp = part.tNext - matAllProp(part.matID).tauImp(part.mode)*log(1-rand());
                    part.relaxCount = part.relaxCount +1;
                end
                
                % resetting starting point
                part.startPoint = part.endPoint;
                
                % finding new direction
                % note that part vel has not been touched yet. Make
                % sure that it is indeed the case.
                
                % finding normal for reflection
                if(dot(part.vel,obj.normal)>0)
                    normalToUse = -obj.normal;
                else
                    normalToUse = obj.normal;
                end
                
                if(rSpec < spec)
                    % specular refelction
                    reflectedVel = Utils.Reflection(part.vel, normalToUse);
                    partDir = reflectedVel/vecnorm(reflectedVel,2,1);
                else
                    % diffusive reflection
                    randDir = Utils.Draw_random('hemisphere');
                    rotMat = Utils.Orient([1;0;0],normalToUse);
                    partDir = rotMat*randDir;
                end
                part.vel = speed*partDir;

            end
            
        end
        
        function nextMode = Find_mode(obj, mode, incMat, isTrans)
            % mode is the initial mode of the particle
            % incMat is the material of incident side 1 or 2 for interface
            % isTrans is true for transmission false for reflection
            if(incMat == 1 && isTrans)
                % transmission from mat1 to mat2
                modalProbs = str2num(obj.TransProbs12{mode});
                indices = str2num(obj.TransIndices12{mode});
            elseif(incMat == 2 && isTrans)
                % transmission from mat2 to mat1
                modalProbs = str2num(obj.TransProbs21{mode});
                indices = str2num(obj.TransIndices21{mode});
            elseif (incMat == 1 && ~isTrans)
                % reflection in mat 1
                modalProbs = str2num(obj.RefProbs12{mode});
                indices = str2num(obj.RefIndices12{mode});
            else
                % reflection in mat2
                modalProbs = str2num(obj.RefProbs21{mode});
                indices = str2num(obj.RefIndices21{mode});
            end
            
            % particle may reflect if there are no triplets that match even
            % momentum for a particular q point. RefProbs will be empty.
            % Temporary fix is that to reflect it with same mode. There
            % will be some error due to resampling of properties
            if(isempty(modalProbs) && ~isTrans)
                nextMode = mode;
                return;
            end
            % drawing random number to compare with cumulative
            % probabilities
            rMode = rand();
            cumulativeSum = cumsum(modalProbs);
            % scaling cumulative sum to match max value to one for choosing
            % correct index
            if(cumulativeSum(end)>0)
                cumulativeSum = cumulativeSum/cumulativeSum(end);
            end
            ii=1;
            while(rMode > cumulativeSum(ii))
                ii = ii + 1;
            end
            % choosing mode corresponding to the chosen probability
            nextMode = indices(ii);
        end
        
        function [success] = Write_flux(obj, directory, identifier)
            fluxFile = [ directory '/FluxAcross' num2str(obj.bndID) '_' identifier '.txt'];
            writematrix(obj.fluxAcross,fluxFile);
            success = true;
        end
    end
end

