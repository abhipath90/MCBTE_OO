classdef DetectorStead < Region
    % DetectorTrans is detector for transient simulations
    %   This has temperature and flux stored in an array of size equal to
    %   the query points in time
    
    properties
        Temp
        Qx
        Qy
        Qz
    end
    
    methods
        function obj = DetectorStead(materials, minX, maxX, minY, maxY, minZ, maxZ, mat)
            % DetectorTrans creates object for transient detector
            %   It is a subclass of Region that contains variables Temp and
            %   3 fluxes as array of lenght equal to time points
            obj@Region(minX, maxX, minY, maxY, minZ, maxZ, mat)
            nmodes = materials(mat).Nmodes;
            obj.Temp = zeros(1,nmodes);
            obj.Qx = zeros(1,nmodes);
            obj.Qy = zeros(1,nmodes);
            obj.Qz = zeros(1,nmodes);
        end
        
        function [obj, Parti, isProblem] = Contribute(obj,Parti,material)
            % Contribute_tans calculates the contribution of particle to
            % this detector's temperature and flux in transient cases
            %   It will check the interaction length between particle and
            %   detector and then calculate its contribution to the
            %   detector.
            isProblem = false;
            len = obj.Len_inside(Parti.startPoint,Parti.endPoint);
            
            if(Parti.matID ~= obj.material && abs(len-0)>2*eps)
%                 partProps = ['Material = ' num2str(Parti.matID) '\nStarting = ' num2str((Parti.startPoint)') ...
%                     '\nEnding = ' num2str((Parti.endPoint)') '\nvelocity = ' num2str((Parti.vel)') ' and \n omega = ' ...
%                     num2str(Parti.omega) '\n'];
                %error(['Material ID of detector and particle do not match\nParticle properties:' partProps]);
                
                %Don't raise error, just kill the particle and count the
                % problematic particles.
                Parti.isAlive = false;
                isProblem = true;
            end
            if(len>0)
                speed = vecnorm(Parti.vel,2,1);
                obj.Temp(1,Parti.mode) = obj.Temp(1,Parti.mode) + Parti.pSign*Parti.eEff*len/material.cvAll/obj.volume/speed;
                obj.Qx(1,Parti.mode) = obj.Qx(1,Parti.mode) + Parti.pSign*Parti.eEff*len*Parti.vel(1)/speed/obj.volume;
                obj.Qy(1,Parti.mode) = obj.Qy(1,Parti.mode) + Parti.pSign*Parti.eEff*len*Parti.vel(2)/speed/obj.volume;
                obj.Qz(1,Parti.mode) = obj.Qz(1,Parti.mode) + Parti.pSign*Parti.eEff*len*Parti.vel(3)/speed/obj.volume;
            end
               
        end
    end
end
