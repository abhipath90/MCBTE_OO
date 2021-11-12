classdef Material
    % This class stores material data
    %   material data including various material data formats, and related
    %   distributions are stored here
    
    properties
        matID
        Nmodes;
        freq
        vel
        tau
        cvMode
        pol
        tauImp
        
        cumulBase
        cumulColl
        cumulVel
        cvAll
    end
    
    methods
        function obj = Material(matInput,ID)
            % Material Construct an instance of this class
            %   Reads raw data and creates Material object which stores
            %   material data and calculates essential distributions
            
            [obj.Nmodes, col] = size(matInput);
            
            obj.matID = ID;
            obj.freq = matInput(:,1);
            obj.vel = matInput(:,2);
            obj.tau = matInput(:,3);
            obj.cvMode = matInput(:,4);
            
            % if polarization is not defined in the material default is
            % same for all branches 1
            if(col>4)
                obj.pol = matInput(:,5);
            else
                obj.pol = ones(size(obj.freq));
            end

            % if impurity scattering times are not defined a very large
            % value is chosen
            if (col>5)
                obj.tauImp = matInput(:,6);
            else
                obj.tauImp = 1e10*ones(obj.Nmodes,1);          % Very large default value
            end
               
            % Calculating distributions
            obj.cumulBase = zeros(obj.Nmodes,1);
            obj.cumulColl = zeros(obj.Nmodes,1);
            obj.cumulVel = zeros(obj.Nmodes,1);
            obj.cumulBase(1) = obj.cvMode(1);
            obj.cumulColl(1) = obj.cvMode(1)/obj.tau(1);
            obj.cumulVel(1) = obj.cvMode(1)*obj.vel(1);
            
            if(obj.Nmodes>1)    % In case gray model is simulated 
                for ii = 2:obj.Nmodes
                    obj.cumulBase(ii) = obj.cumulBase(ii-1) + obj.cvMode(ii);
                    obj.cumulColl(ii) = obj.cumulColl(ii-1) + obj.cvMode(ii)/obj.tau(ii);
                    obj.cumulVel(ii) = obj.cumulVel(ii-1) + obj.cvMode(ii)*obj.vel(ii);
                end
            end
            obj.cvAll = obj.cumulBase(end);
        end
        

    end
end

