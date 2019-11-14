classdef PH_Element < matlab.mixin.Copyable
    properties (Access = public) 
        type
        ports
        attributes
    end
    
    methods (Access = public)
        function obj = PH_Element(type, ports, varargin)
            obj.type = type;
            obj.ports = ports;
            n_attr = floor(nargin - 2)/2;
            obj.attributes = cell(n_attr, 2);
            for a=1:n_attr
                obj.attributes{a, 1} = varargin{2*a-1};
                obj.attributes{a, 2} = varargin{2*a};
            end
        end
    end
end