classdef (Abstract) PH_System < matlab.mixin.Copyable
    %PORTHAMILTONIANSYSTEM An absctract superclass for a port-Hamiltonian
    %System
    %   Defines properties and methods for a very general port-Hamiltonian system 
    
    properties (Access = public)
        % name of the port-Hamiltonian system    
        name = 'port-Hamiltonian system'
    end
    
    properties (SetAccess = protected, GetAccess = public) 
        n           % System order
        n_u         % Number of inputs/external/boundary ports
        n_c         % Number of coupling constraints
        
        % Description of the system's ports
        n_ports     % Number of ports
        ports       % Contains information on the system's ports
        
        % Description of the system's elements
        n_elements  % Number of elements
        elements    % System elements
        
        % System nodes 
        n_nodes     % Number of nodes
        nodes       % System nodes
    end
    
    methods (Access = public, Abstract)
        y       = getSystemOutput(obj, x, u)
      
        add(obj, system)            % Add another element or subsystem
    end
    
    methods (Access = public)
        % Constructor
        function obj = PH_System(name)
            if nargin > 0
                if ~isa(name, 'char')
                    error('PHSystem: name must be given as a string');
                end
                obj.name = name;
            end
            
            obj.n_ports = 0;
            obj.ports = cell(0); 
            obj.n_elements = 0;
            obj.elements = cell(0);
            obj.n_nodes = 0;
            obj.nodes = cell(0);
        end
        
        % System properties
        function ret = isConstrained(obj)
            ret = obj.n_c > 0;
        end
    end
end

