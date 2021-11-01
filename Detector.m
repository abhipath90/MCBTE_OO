classdef Detector
    %Detector stores all the detectors present in the simulation
    %   The detectors are setup based on the simulation type and domain
    %   definition and thermodynamic variables are calculated
    
    properties
        detectors       % array of detector objects
        noDetect        % number of detectors
    end
    
    methods
        function obj = Detector(SimObj,measureReg,measureTime,matProps)
            % Detector is Construct an instance of this class
            %   Detailed explanation goes here
            [obj.noDetect,~] = size(measureReg);
            obj.detectors = [];
            if(strcmpi(SimObj.typeSim,'Steady'))
                % Generate steady state detectors
                for ii=1:obj.noDetect
                    minx = measureReg(ii,1); maxx = measureReg(ii,2);
                    miny = measureReg(ii,3); maxy = measureReg(ii,4);
                    minz = measureReg(ii,5); maxz = measureReg(ii,6);
                    mat = measureReg(ii,7);
                    obj.detectors = [obj.detectors; DetectorStead(matProps,minx,maxx,miny,maxy,minz,maxz,mat)];
                end
            else
                % Generate transient detectors
                for ii=1:obj.noDetect
                    minx = measureReg(ii,1); maxx = measureReg(ii,2);
                    miny = measureReg(ii,3); maxy = measureReg(ii,4);
                    minz = measureReg(ii,5); maxz = measureReg(ii,6);
                    mat = measureReg(ii,7);
                    obj.detectors = [obj.detectors; DetectorTrans(minx,maxx,miny,maxy,minz,maxz,mat,measureTime)];
                end
            end
            
        end
        
        function obj = Record_contributions(obj,part,matProp)
            % Record_contributions records particle contribution to proper
            % detectors
            
            % The detectors have methods with exactly same names, based on
            % their existance in the simulation, correct method will
            % automatically be called. No need to identify type of detector
            % here
            
            for ii=1:obj.noDetect
                obj.detectors(ii) = obj.detectors(ii).Contribute(part,matProp(part.matID));
            end
                        
        end
        
        function [success] = Write_output(obj, directory, identifier)
            success = false;
            tempFile = [ directory '/Temp' identifier '.txt'];
            qxFile = [ directory '/Qx' identifier '.txt'];
            qyFile = [ directory '/Qy' identifier '.txt'];
            qzFile = [ directory '/Qz' identifier '.txt'];
            
            tempmat = zeros(obj.noDetect,length(obj.detectors(1).Temp));
            qxmat = zeros(size(tempmat));
            qymat = zeros(size(tempmat));
            qzmat = zeros(size(tempmat));
            for ii=1:obj.noDetect
                tempmat(ii,:) = obj.detectors(ii).Temp;
                qxmat(ii,:) = obj.detectors(ii).Qx;
                qymat(ii,:) = obj.detectors(ii).Qy;
                qzmat(ii,:) = obj.detectors(ii).Qz;
            end
            writematrix(tempmat, tempFile);
            writematrix(qxmat, qxFile);
            writematrix(qymat, qyFile);
            writematrix(qzmat, qzFile);
            
            success = true;
        end
    end
end

