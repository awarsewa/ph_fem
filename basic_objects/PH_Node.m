classdef PH_Node < matlab.mixin.Copyable
    properties (Access = public) 
        domain      % Physical domain of this node (e.g. mechanical/electrial/...)
        ports       % Numbers of the ports connected to this node
        elements    % Numbers of the elments connected to this node
    end
    
    % TODO: error checking
    methods (Access = public)
        function obj = PH_Node(domain, ports, elements)
            obj.domain = domain;
            if nargin > 1
                obj.ports = ports;
            else
                obj.ports = [];
            end
            if nargin > 2
                obj.elements = elements;
            else
                obj.elements = [];
            end
        end
    end
end