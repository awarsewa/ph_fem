classdef PH_Port < matlab.mixin.Copyable
    properties (Access = public) 
        scope           % storage/external
        domain          % Physical domain of this port (e.g. mechanical/electrial/...)
        type            % Port type (e.g. mechanical -> translation/rotation/translatory strain/rotational strain)
        nodes           % Nodes to which the port is connected 
    end
    
    methods (Access = public)
        function obj = PH_Port(scope, domain, type, nodes)
            if ~ischar(scope) || ~(strcmp(scope, 'storage') || strcmp(scope, 'external') || strcmp(scope, 'boundary'))
                error('domain must be either ''storage'', ''external'' or ''boundary''.');  
            end
            
            if ~ischar(domain) || ~(strcmp(domain, 'mechanical') || strcmp(domain, 'electrical') || strcmp(domain, 'hydraulic'))
                error('domain must be either ''electrical'', ''mechanical'' or ''hydraulic''.');
            end
            
            
            if ~ischar(type)
                error('type must be supplied as a character array');
            end
            
            if ~iscolumn(nodes)
                error('port nodes must be supplied as a column vector');
            end
            
            obj.domain = domain;
            obj.scope = scope;
            obj.type = type;
            obj.nodes = nodes;
        end
    end
end