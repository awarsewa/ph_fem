classdef PH_MechanicalPort < PH_Port
    properties (Access = public) 
        orientation
    end
    
    methods (Access = public)
        function obj = PH_MechanicalPort(scope, type, nodes, orientation)
            % Type relates to input type for a mechanical port
            if ~ischar(type) || ~(  strcmp(type, 'force') || ...
                                    strcmp(type, 'torque') || ...
                                    strcmp(type, 'velocity') || ...
                                    strcmp(type, 'angular velocity'))        
                error('type must be supplied as one of ''force'', ''torque'', ''velocity'' or ''angular velocity''.');
            end
            
            obj = obj@PH_Port(scope, 'mechanical', type, nodes);
            
            if ~ismatrix(orientation) || length(nodes) ~= size(orientation, 2) || size(orientation, 1) ~= 3
                error('orientation must be supplied as a 3xn matrix, where n is the number of nodes the port is connected to');
            end
        end
    end
end