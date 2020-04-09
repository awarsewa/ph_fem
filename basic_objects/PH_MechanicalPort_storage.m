classdef PH_MechanicalPort_storage < PH_MechanicalPort & PH_Port_storage
    methods (Access = public)
        function obj = PH_MechanicalPort_storage(type, nodes, orientation, effort, flow)        
            obj@PH_MechanicalPort('storage', type, nodes, orientation);
            obj@PH_Port_storage('mechanical', type, nodes, effort, flow);
            obj.orientation = orientaion;
            obj.effort = effort;
            obj.flow = flow;
        end
    end
end