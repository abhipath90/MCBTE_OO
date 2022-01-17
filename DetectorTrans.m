classdef DetectorTrans < Region
    % DetectorTrans is detector for transient simulations
    %   This has temperature and flux stored in an array of size equal to
    %   the query points in time
    
    properties
        Temp
        Qx
        Qy
        Qz
        timePoints
    end
    
    methods
        function obj = DetectorTrans(minX, maxX, minY, maxY, minZ, maxZ, mat,timePoints)
            % DetectorTrans creates object for transient detector
            %   It is a subclass of Region that contains variables Temp and
            %   3 fluxes as array of lenght equal to time points
            obj@Region(minX, maxX, minY, maxY, minZ, maxZ, mat)
            obj.timePoints = timePoints;
            times = length(timePoints);
            obj.Temp = zeros(1,times);
            obj.Qx = zeros(1,times);
            obj.Qy = zeros(1,times);
            obj.Qz = zeros(1,times);
        end
        
        function [obj,Parti,isProblem] = Contribute(obj,Parti,material)
            % Contribute_tans calculates the contribution of particle to
            % this detector's temperature and flux in transient cases
            %   It will check if there is an interaction between particle
            %   flight and the detector at the time of interest and then
            %   calculate its contribution
            
            isProblem = false;
            % find timeIdexes
            timeIndex = find(obj.timePoints>Parti.t0 & obj.timePoints<Parti.tNext);
            
            for ii=1:length(timeIndex)
                timeStamp = obj.timePoints(timeIndex(ii));
                
                Position = Parti.startPoint + Parti.vel*(timeStamp - Parti.t0);
                isInside = obj.Check_inside(Position);
                
                if(isInside)
                    if(Parti.matID ~= obj.material)
                        
                        Parti.isAlive = false;
                        isProblem = true;
                        %error('Material ID of detector and particle do not match');
                    end
                    %matID = Parti.matID;
                    obj.Temp(1,timeIndex(ii)) = obj.Temp(1,timeIndex(ii)) + Parti.pSign*Parti.eEff/material.cvAll/obj.volume;
                    obj.Qx(1,timeIndex(ii)) = obj.Qx(1,timeIndex(ii)) + Parti.pSign*Parti.eEff*Parti.vel(1)/obj.volume;
                    obj.Qy(1,timeIndex(ii)) = obj.Qy(1,timeIndex(ii)) + Parti.pSign*Parti.eEff*Parti.vel(2)/obj.volume;
                    obj.Qz(1,timeIndex(ii)) = obj.Qz(1,timeIndex(ii)) + Parti.pSign*Parti.eEff*Parti.vel(3)/obj.volume;
                end
            end
               
        end
    end
end

