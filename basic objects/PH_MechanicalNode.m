classdef PH_MechanicalNode < PH_Node
    % Mechanical ports have an additional orientation in x/y/z coordinates
    properties (Access = public) 
        location    % Location of this port in x/y/z coordinates
        lockedDOFs  % DOFs that are locked/fixed
    end
    methods (Access = public)
        function obj = PH_MechanicalNode(ports, location, lockedDOFs)
            obj = obj@PH_Node('mechanical', ports);
            
            if isempty(location) || ~isrow(location) || length(location)  ~= 3
                error('location must be supplied as a 1x3 vector')
            end
            obj.location = location;
            
            if nargin > 2
                if isempty(lockedDOFs) || ~isrow(lockedDOFs) || ~(length(lockedDOFs) == 6)
                    error('lockedDOFs must be supplied as a 1x6 vector');
                end
                obj.lockedDOFs = lockedDOFs;
            else
                obj.lockedDOFs = zeros(1, 6);
            end
        end
    end
end