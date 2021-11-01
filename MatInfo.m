classdef MatInfo
    % MatInfo class stores the information about the material properties
    %   Here we store information about the characteristics of the
    %   materials. This includes their number, type, name etc.
    
    properties
        matNum          % number of materials in the simulation
        matName         % Array of material names
        matProp          % Material object arrary
    end
    
    methods
        function obj = MatInfo(RawData)
            %UNTITLED15 Construct an instance of this class
            %   Detailed explanation goes here
            [obj.matNum,~]= size(RawData.matData);
            obj.matName = [];
            obj.matProp = [];
            for ii=1:obj.matNum
                obj.matName = [obj.matName; ['material' num2str(ii)] ];
                obj.matProp = [obj.matProp; Material(RawData.matData{ii},ii)];
            end
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

