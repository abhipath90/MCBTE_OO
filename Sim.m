classdef Sim
    % This contains information about the simulation parameters
    %  
    
    properties
        N           % Number of particles
        Teq         % Equilibirum temperaure
        maxScatter  % Maximum scattering events allowed
        tMax        % Maximum simulation time
        typeSim     % Type of simulations Steady/transient
        vol         % Volume of the domain (currently for gradient source)
    end
    
    methods
        function obj = Sim(RawData)
            % RawData is an InputData object that contains all the
            % information that is read from the files
            obj.N = RawData.param(1);
            obj.Teq = RawData.param(4);
            obj.maxScatter = RawData.param(2);
            
            if (isempty(RawData.measureTime))
                maxTime = [];
                simType = 'Steady';
            else
                maxTime = RawData.measureTime(end);
                simType = 'Transient';
            end
            
            obj.tMax = maxTime;
            obj.typeSim = simType;
            obj.vol = RawData.param(3);
        end
        
    end
end

