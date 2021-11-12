classdef Region
    % Region class describes a geometric region in the domain
    %   Curretnly it is only supporting coordinate extents as it is
    %   primarily used to create detectors
    
    properties
        xmax
        xmin
        ymax
        ymin
        zmax
        zmin
        volume
        material    % material ID
    end
    
    methods
        function obj = Region(minX, maxX, minY, maxY, minZ, maxZ, mat)
            %Region Construct an instance of this class
            %   This will assign values to the properties defined
            obj.xmin = minX;
            obj.xmax = maxX;
            obj.ymin = minY;
            obj.ymax = maxY;
            obj.zmin = minZ;
            obj.zmax = maxZ;
            obj.volume = (maxX-minX)*(maxY-minY)*(maxZ-minZ);
            obj.material = mat;
            
        end
        
        function isInside= Check_inside(obj,Position)
            % Check_inside checks if particle belongs inside a region at a
            % time given
            %   This is useful in finding particle contribution to flux and
            %   temperature in transient case
            isInside = false;
            Xc = Position(1); Yc = Position(2); Zc = Position(3);
            X_prod = (Xc-obj.xmin)*(Xc-obj.xmax);
            Y_prod = (Yc-obj.ymin)*(Yc-obj.ymax);
            Z_prod = (Zc-obj.zmin)*(Zc-obj.zmax);
            
            if(X_prod<0 && Y_prod<0 && Z_prod<0)
                isInside = true;
            end
        end
        
        function len = Len_inside(obj,Start,End)
            
        % Note that tmin and tmax are not parametric fractions
        % Their value is not in [0,1] for proper interaction
        % it is in [0,lenSeg] for proper interaction
            
            len = 0;
            lenSeg = vecnorm(End-Start,2,1);

            cosine_vec = (End-Start)/lenSeg; % direction cosines
            x0 = Start(1); y0 = Start(2); z0 = Start(3);
            lx = cosine_vec(1); ly = cosine_vec(2); lz = cosine_vec(3);
            % Finding if the detector is in the path of the flight
            
            % Finding parametric fractions when the particle will cross X
            % boundaries
            if(abs(lx-0)<2*eps)
                % particle is running perpendicular to x axis, its x
                % coordinate doesn't change
                if((obj.xmin-x0)*(obj.xmax-x0)<0)
                    % particle started inside the bounds, will remain there
                    txmin = 0;
                    txmax = lenSeg;
                else
                    % particle started outside the x bounds and will
                    % remain there. No interaction
                    return;
                end
            elseif (lx>0)
                % x coordinate increases with time
                txmin = (obj.xmin - x0)/lx;
                txmax = (obj.xmax - x0)/lx;
            else
                % x coordinate decreases with time
                txmin = (obj.xmax - x0)/lx;
                txmax = (obj.xmin - x0)/lx;
                
            end
            
            tmin = txmin;
            tmax = txmax;
            
            if(abs(ly-0)<2*eps)
                % particle is running perpendicular to the y axis. its y
                % cooridnate doesn't change
                if((obj.ymin-y0)*(obj.ymax-y0)<0)
                    % Particle start within the bounds and will remain
                    % inside those bounds
                    tymin = 0;
                    tymax = lenSeg;
                else
                    % particle starts outside of the bounds and will remain
                    % outside the bounds. No interaction
                    return;
                end
            elseif (ly>0)
                % y coordinate increase with time
                tymin = (obj.ymin - y0)/ly;
                tymax = (obj.ymax - y0)/ly;
            else
                % y coordinate decrease with time
                tymin = (obj.ymax - y0)/ly;
                tymax = (obj.ymin - y0)/ly;
                
            end
            
            if((tmin>tymax) || (tymin > tmax))
                % fractions during which particle stays in x bounds ane
                % exclusive of fractions for which particle stays in y
                % bounds. No interaction
                return;
            end
            
            % adjusting the bounds 
            if (tymin > tmin)
                tmin = tymin;
            end
            if (tymax < tmax)
                tmax = tymax;
            end
            
            if(abs(lz-0)<2*eps)
                % particle is running perpendicular to the z axis. its z
                % coordinate doesn't change
                if((obj.zmin - z0)*(obj.zmax-z0)<0)
                    % particle starts within the bounds and will remain
                    % within the bounds
                    tzmin = 0;
                    tzmax = lenSeg;
                else
                    % particle starts outside the bounds and will remain
                    % outside the bounds
                    return;
                end
            elseif (lz>0)
                % z coordinate increase with time
                tzmin = (obj.zmin - z0)/lz;
                tzmax = (obj.zmax - z0)/lz;
            else
                % z coordinate decreases with time
                tzmin = (obj.zmax - z0)/lz;
                tzmax = (obj.zmin - z0)/lz;
            end
            
            if((tmin > tzmax) || (tmax < tzmin))
                % fractions during which particle stays in z bounds are
                % exclusive of the fractions during which particle stays in
                % the x-y bounds
                return;
            end
            
            % adjusting the bounds
            if(tzmin > tmin)
                tmin = tzmin;
            end
            if(tzmax < tmax)
                tmax = tzmax;
            end
            
            % If we reach here, the detector is in the path of the particle
            % flight. A real interaction will happen only if the parameters
            % calculated are in [0,lenSeg] after normalizing with the lenght of
            % the flight
            
            % Calculating length of the segment
            if (tmin > lenSeg || tmax <0)
                % The flight is completely outside the bounds
                return;
            end
            
            if((tmin<=(0+2*eps)) && (tmax >=lenSeg))
                % The flight is completely inside the bounds
                len = lenSeg;
                return;
            end
            
            if((tmin>0) && (tmax<lenSeg))
                % flight enters and then exits the bounds
                len = tmax - tmin;
                return;
            elseif ( (tmin>0) && (tmax >=lenSeg))
                % flight enters the bounds but doesn't go out
                len = lenSeg - tmin;
                return;
            elseif ((tmin<=0) && (tmax < lenSeg))
                % flight starts inside and then goes out of bounds
                len = tmax;
                return;
            else
                disp(['start ' num2str(Start) 'end ' num2str(End)])
                error('Error while calculating flight-Region interaction.');
            end
        end
        
        function randPoint = Random_point(obj)
            xCord = obj.xmin + rand()*(obj.xmax-obj.xmin);
            yCord = obj.ymin + rand()*(obj.ymax-obj.ymin);
            zCord = obj.zmin + rand()*(obj.zmax-obj.zmin);
            randPoint = [xCord; yCord; zCord];
        end
        
        
        
    end
end

