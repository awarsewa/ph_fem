classdef PH_Port_storage < PH_Port
    properties (Access = public) 
        effort  % number of the port's effort
        flow    % number of the port's flow
    end
    
    methods (Access = public)
        function obj = PH_Port_storage(domain, type, nodes, effort, flow)
            if ~isscalar(effort) || ~isscalar(flow)
                error('flow and effort must be supplied as a row vectors'); 
            end
            
            % Instantiate superclass
            obj@PH_Port('storage', domain, type, nodes);
            
            obj.effort = effort;
            obj.flow = flow;
        end
    end
end