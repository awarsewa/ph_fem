classdef PH_Port_boundary < PH_Port
    properties (Access = public) 
        element     % Identifies the element this port is associated with
        IOPair      % Identifies the input/output pair associated with this boundary port
        inputName   % Name of the input (e.g. n1_Fx for a force in x-direction acting on node 1)
        outputName  % Name of the output (e.g. n1_vz for the corresponding velocity
    end
    
    methods (Access = public)
        function obj = PH_Port_boundary(domain, type, nodes, element, IOPair, inputName, outputName)
            if ~isscalar(element)
                error('element nr. must be supplied as a scalar');
            end
            if ~isscalar(IOPair)
                error('input/output nr. must be supplied as a scalar'); 
            end
            if ~ischar(inputName) || ~ischar(outputName)
                error('input/output names must be supplied as character arrays');
            end
            
            obj@PH_Port('boundary', domain, type, nodes); 
            
            obj.element = element;
            obj.IOPair = IOPair;
            obj.inputName = inputName;
            obj.outputName = outputName;
        end
    end
end