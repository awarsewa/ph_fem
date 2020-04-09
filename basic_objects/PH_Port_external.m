classdef PH_Port_external < PH_Port
    properties (Access = public) 
        IOPair      % Identifies the input/output pair associated with this boundary port
        inputName   % Name of the input (e.g. n1_Fx for a force in x-direction acting on node 1)
        outputName  % Name of the output (e.g. n1_vz for the corresponding velocity
    end
    
    methods (Access = public)
        function obj = PH_Port_external(domain, type, nodes, IOPair, inputName, outputName)
            if ~isscalar(IOPair)
                error('input/output nr. must be supplied as a scalar'); 
            end
            if ~ischar(inputName) || ~ischar(outputName)
                error('input/output names must be supplied as character arrays');
            end
            
            obj@PH_Port('external', domain, type, nodes); 
            
            obj.IOPair = IOPair;
        end
    end
end