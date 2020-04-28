classdef PH_MechanicalPort_boundary < PH_MechanicalPort & PH_Port_boundary
    methods (Access = public)
        function obj = PH_MechanicalPort_boundary(type, nodes, orientation, element, IOPair, inputName, outputName)        
            obj@PH_MechanicalPort('boundary', type, nodes, orientation);
            obj@PH_Port_boundary('mechanical', type, nodes, element, IOPair, inputName, outputName);
            obj.element = element;
            obj.inputName = inputName;
            obj.outputName = outputName;
            obj.IOPair = IOPair;
            obj.orientation = orientation;
        end
    end
end