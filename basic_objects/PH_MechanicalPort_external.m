classdef PH_MechanicalPort_external < PH_MechanicalPort & PH_Port_external
    methods (Access = public)
        function obj = PH_MechanicalPort_external(type, nodes, orientation, IOPair, inputName, outputName)        
            obj@PH_MechanicalPort('external', type, nodes, orientation);
            obj@PH_Port_external('mechanical', type, nodes, IOPair, inputName, outputName);
            obj.inputName = inputName;
            obj.outputName = outputName;
            obj.IOPair = IOPair;
            obj.orientation = orientation;
        end
    end
end