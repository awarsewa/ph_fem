classdef PH_Element < matlab.mixin.Copyable
    properties (Access = public) 
        type        % element type (e.g. spring)
        nodes       % nodes the element is connected to
        ports       % numbers of the ports associated with the element
        attributes  % additional attributes of the element
    end
    
    % TODO: error checking
    methods (Access = public)
        function obj = PH_Element(type, nodes, ports, varargin)
            obj.type = type;
            obj.nodes = nodes;
            obj.ports = ports;
            n_attr = floor(nargin - 3)/2;
            obj.attributes = cell(n_attr, 2);
            for a=1:n_attr
                obj.attributes{a, 1} = varargin{2*a-1};
                obj.attributes{a, 2} = varargin{2*a};
            end
        end
    end
end