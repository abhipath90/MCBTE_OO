classdef Geometry
    % Geometry class contains all details about the geometry of the domain
    %   Number of boundaries of each type, their object array. There will
    %   be methods as well to help particle move forward in time. The first
    %   one will be the one to find the closest intraction 
    
    properties
        numBnd          % Number of boundaries in the domain
        numIso          % Number of isothermal boundaries
        numAdi          % Number of adiabatic boundaries
        numPeri         % Number of periodic boundaries
        numInter        % Number of interfaces
        Iso             % Array of isothermal boundaries
        Adi             % Array of adiabatic boundaries
        Peri            % Array of periodic boundaries
        Inter           % Array of interface boundaries
    end
    
    methods
        function obj = Geometry(RawData,matObj)
            % Geometry constructor from raw data
            %   This method reads raw data and initializes all the boundary
            %   objects to finishe the definition of the geometry
            
            % finding number of all the boundaries defining the domain
            bndData = RawData.bnd;
            [obj.numBnd, ~] = size(bndData);
            
            % initializing isothermal boundaries
            obj.Iso = [];
            obj.numIso = 0;
            % Finding isothermal boundaries
            indexIso = find(RawData.bndProp(:,2)==1);
                        
            if( ~isempty(indexIso))
                obj.numIso = length(indexIso);
                for ii=1:obj.numIso
                    thisBnd = RawData.bnd(indexIso(ii),:);
                    thisBndProp = RawData.bndProp(indexIso(ii),:);
                    point1 = thisBnd(1,1:3);
                    point2 = thisBnd(1,4:6);
                    point3 = thisBnd(1,7:9);
                    point4 = thisBnd(1,10:12);
                    mat = thisBndProp(1,3);
                    temp = thisBndProp(1,4);
                    obj.Iso = [obj.Iso ; IsoBnd(indexIso(ii),mat,point1,point2,point3,point4,temp)];
                end
            end
            
            % initializing adiabatic boundaries
            obj.Adi = [];
            obj.numAdi = 0;
            % finding adiabatic boundaries
            indexAdi = find(RawData.bndProp(:,2)==2);
            
            if(~isempty(indexAdi))
                obj.numAdi = length(indexAdi);
                for ii=1:obj.numAdi
                    thisBnd = RawData.bnd(indexAdi(ii),:);
                    thisBndProp = RawData.bndProp(indexAdi(ii),:);
                    point1 = thisBnd(1,1:3);
                    point2 = thisBnd(1,4:6);
                    point3 = thisBnd(1,7:9);
                    point4 = thisBnd(1,10:12);
                    mat = thisBndProp(1,3);
                    spec = thisBndProp(1,4);
                    obj.Adi = [obj.Adi; AdiaBnd(indexAdi(ii),mat,point1,point2,point3,point4,spec)];
                end
            end
            
            % initializing periodic boundaries
            obj.Peri = [];
            obj.numPeri = 0;
            % Finding periodic boundaries
            indexPeri = find(RawData.bndProp(:,2)==3);
            bndPairs = RawData.periPair;
            if(~isempty(indexPeri))
                obj.numPeri = length(indexPeri);
                for ii=1:obj.numPeri
                    thisBnd = RawData.bnd(indexPeri(ii),:);
                    thisBndProp = RawData.bndProp(indexPeri(ii),:);
                    point1 = thisBnd(1,1:3);
                    point2 = thisBnd(1,4:6);
                    point3 = thisBnd(1,7:9);
                    point4 = thisBnd(1,10:12);
                    mat = thisBndProp(1,3);
                    transVec = (thisBndProp(1,4:6))'; % transposing for column vector
                    pairBnd = bndPairs(bndPairs(:,1) == indexPeri(ii),2);
                    obj.Peri = [obj.Peri; PeriBnd(indexPeri(ii),mat,point1,point2,point3,point4,transVec,pairBnd)];
                end
            end
            
            % initializing interface boundaries
            obj.Inter = [];
            obj.numInter = 0;
            % Finding interface boundaries in the simulation
            indexInter = find(RawData.bndProp(:,2)==5);
            
            if(~isempty(indexInter))
                obj.numInter = length(indexInter);
                if(obj.numInter ~= RawData.interface{1})
                    error('number of interface definitions do not match');
                end
                % Ordered array of all the interfaces listed in interface
                % data
                bndArray = RawData.interface{2};
                % Corresponding interface characteristic data filenames
                dataPoints = RawData.interface{3};
                for ii=1:obj.numInter
                    % Going through the list of all the boundaries in
                    % simulation
                    thisBnd = RawData.bnd(indexInter(ii),:);
                    thisBndProp = RawData.bndProp(indexInter(ii),:);
                    point1 = thisBnd(1,1:3);
                    point2 = thisBnd(1,4:6);
                    point3 = thisBnd(1,7:9);
                    point4 = thisBnd(1,10:12);
                    mat = thisBndProp(1,3);
                    mat2 = thisBndProp(1,4);
                    % To find where the data for this boundary is located
                    % in the interface data filelist
                    indexData = bndArray == indexInter(ii);
                    % List of files that define interface characteristics.
                    type = dataPoints{indexData};
                    obj.Inter = [obj.Inter; InterBnd(indexInter(ii),mat,point1,point2,point3,point4,mat2,type,matObj)];
                end
            end
        end
        
        function [time,point,bndID]= Geometric_interaction(obj,part)
            % Geometric_interaction calculates closest interaction in
            % future
            %   It circles through all the boundaries defined in the
            %   geometry and finds out the boundary that is closest and at
            %   what time the interaction will happen
            frac = zeros(obj.numBnd,1); point = zeros(3,obj.numBnd);
            bndID = zeros(obj.numBnd,1);
            for ii=1:obj.numBnd
                if(ii<=obj.numIso)
                    % checking interaction with isothermal walls
                    [frac(ii),point(:,ii),bndID(ii)] = obj.Iso(ii).interaction(part);
                elseif(ii<=(obj.numIso+obj.numAdi))
                    % checking interaction with adiabatic walls
                    subtract = obj.numIso;
                    [frac(ii),point(:,ii),bndID(ii)] = obj.Adi(ii-subtract).interaction(part);
                elseif(ii<=(obj.numIso+obj.numAdi+obj.numPeri))
                    % checking interaction with periodic walls
                    subtract = obj.numIso+obj.numAdi;
                    [frac(ii),point(:,ii),bndID(ii)] = obj.Peri(ii-subtract).interaction(part);
                else
                    subtract = obj.numIso + obj.numAdi + obj.numPeri;
                    [frac(ii),point(:,ii),bndID(ii)] = obj.Inter(ii-subtract).interaction(part);
                end
            end
                
            % finding proper intreactions
            index = find(frac>0);
            
            if(isempty(index))
                % no interaction
                time = []; point = []; bndID = [];
                return;
            else
                frac = frac(index,:);
                point = point(:,index); % point contains points in columns
                bndID = bndID(index,:);
                
                [~,minIndex] = min(frac);
                
                frac = frac(minIndex);
                point = point(:,minIndex);
                bndID = bndID(minIndex);
                
                time = part.t0 + frac*(part.tNext-part.t0);
            end
        end
        
        
        function Draw(obj)
            % this method is to draw the whole geometry with boundary
            % normals
            
            figure()
            
            for ii=1:obj.numBnd
                
                if(ii<=obj.numIso)
                    bnd = obj.Iso(ii);
                elseif(ii<=(obj.numIso+obj.numAdi))
                    subtract = obj.numIso;
                    bnd = obj.Adi(ii-subtract);
                elseif(ii<=(obj.numIso+obj.numAdi+obj.numPeri))
                    subtract = obj.numIso+obj.numAdi;
                    bnd = obj.Peri(ii-subtract);
                else
                    subtract = obj.numIso + obj.numAdi + obj.numPeri;
                    bnd = obj.Inter(ii-subtract);
                end
                xs = [bnd.vertex1(1);bnd.vertex2(1);bnd.vertex3(1);bnd.vertex4(1)];
                ys = [bnd.vertex1(2);bnd.vertex2(2);bnd.vertex3(2);bnd.vertex4(2)];
                zs = [bnd.vertex1(3);bnd.vertex2(3);bnd.vertex3(3);bnd.vertex4(3)];
                midpoint = (bnd.vertex1+bnd.vertex3)/2;
                diagLen = vecnorm(bnd.vertex1-bnd.vertex3,2,1);
                diff = bnd.normal*0.3*diagLen;
                fill3(xs,ys,zs,'g','FaceAlpha',0.5);
                if(ii==1)
                    hold on;
                end
                quiver3(midpoint(1),midpoint(2),midpoint(3),diff(1),diff(2),diff(3));
            end
            view(30,30);
            xlabel('X'); ylabel('Y'); zlabel('Z');
        end
        
    end
end

