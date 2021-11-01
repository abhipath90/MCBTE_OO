classdef Utils
    % Utilities contains some useful methods used by all classes
    %   Detailed explanation goes here
    
    properties
        
    end
    
    methods (Static)
        
        function outVec = Reflection(inVec,normal)
            % this method implements householder transform
            outVec = inVec - 2*normal*(normal'*inVec)/(normal'*normal);
        end
        
        function R = Orient(ref, dest)
            
            % This code generates the rotation matrix to rotate a vector ref to vector dest
            % This is used many times to orient the velocity vector appropriately when
            % emitted or diffusively reflected from a plane
            % WORKS FOR ANY DIMENSION
            
            % This procedure only works when ref~=-(dest);
            
            if(~all(ref==-dest))
                u = ref/vecnorm(ref,2,1);
                v = dest/vecnorm(dest,2,1);
                N = length(u);
                S = Utils.Reflection(eye(N),v+u); % SU=-v; Sv=-u;
                R = Utils.Reflection(S,v);        % v=Su;
            else
                R = -eye(length(ref));
            end
        end
        
        function index = Select_mode(cumul,Nmodes)
            %SELECT_MODE - This function randomly draws an index within a cumulative
            %distribution. The distribution of the indexes follows the probability
            %density function from which the cumulative distribution is generated.
            %   The cumulative distribution cumul must be a 1D array with increasing
            %   values. Nmodes is the length of the array.
            
            R = rand();
            i1 = 0;
            i3 = Nmodes;
            i2 = floor((i1+i3)/2);
            
            while (i3-i1>1)
                
                if R<cumul(i2)/cumul(Nmodes)
                    i3  = i2;
                    i2 = floor((i1+i3)/2);
                else
                    i1 = i2;
                    i2 = floor((i1+i3)/2);
                end
            end
            index = i3;
            
        end
        
        function randDir = Draw_random(type)
            if(strcmpi(type,'hemisphere'))
                % drawing random direction on hemisphere
                theta = 2*pi*rand();    % azimuthal angle
                u = rand();             % for polar angle
                cosPhi = sqrt(u);
                sinPhi = sqrt(1-u);
                randDir = [cosPhi; sinPhi*cos(theta); sinPhi*sin(theta)];
                
            elseif(strcmpi(type,'sphere'))
                % drawing random direction on sphere
                theta = 2*pi()*rand();  % azimuthal angle
                u = rand();             % for sampling polar angle
                cosPhi = 1-2*u;
                sinPhi = sqrt(1-cosPhi^2);
                randDir = [cosPhi; sinPhi*cos(theta); sinPhi*sin(theta)];
            else
                error('unexpected input for choosing random direction');
            end
        end

    end
end

