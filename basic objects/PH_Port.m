classdef PH_Port < matlab.mixin.Copyable
    properties (Access = public) 
        domain      % Physical domain of this port (e.g. mechanical/electrial/...)
        scope       % External/internal/boundary port
        type        % Port type (e.g. mechanical -> force/velocity
        IOPair      % input/output number associated with this port
        inputType   % For a force/velocity port, is the input force or velocity? 
        outputType  % Same with respect to the output
        node        % Node to which the port is connected 
    end
    
    methods (Access = public)
        function obj = PH_Port(domain, scope, type, IOPair, inputType, outputType, node)
            obj.domain = domain;
            obj.scope = scope;
            obj.inputType = inputType;
            obj.outputType = outputType;
            obj.type = type;
            obj.IOPair = IOPair;
            obj.node = node;
        end
    end
end