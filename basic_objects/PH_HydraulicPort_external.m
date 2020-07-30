classdef PH_HydraulicPort_external < PH_Port_external  
    methods (Access = public)
        function obj = PH_HydraulicPort_external(type, IOPair, inputName, outputName)
            % Type relates to input type for a mechanical port
            if ~ischar(type) || ~strcmp(type, 'valve')        
                error('type must be ''valve'' atm.');
            end
            
            obj@PH_Port_external('hydraulic', type, 0, IOPair, inputName, outputName);
            obj.inputName = inputName;
            obj.outputName = outputName;
            obj.IOPair = IOPair;
        end
    end
end