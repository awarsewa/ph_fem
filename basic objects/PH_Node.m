classdef PH_Node < matlab.mixin.Copyable
    properties (Access = public) 
        domain      % Physical domain of this node (e.g. mechanical/electrial/...)
        ports
    end
    
    methods (Access = public)
        function obj = PH_Node(domain, ports)
            obj.domain = domain;
            if nargin > 1
                obj.ports = ports;
            else
                obj.ports = [];
            end
        end
    end
end