classdef (Abstract) PH_FEM_Beam_element < PH_FEM_element 
    methods(Access = public)
        function obj = PH_FEM_Beam_element(n, J, Q, varargin)
            obj = obj@PH_FEM_element(n, J, Q, varargin{:});
        end
        
        function setPosition(obj, nodePositions)
            obj.setPosition@PH_FEM_element(nodePositions);
        end
        
        function setAttributes(obj, attribs)
            if isempty(attribs{1})
                return
            end
            
            % First attribute of a beam element is the orientation of its
            % bending axis
            axis = attribs{1};
            if ~iscolumn(axis)
                axis = axis';
            end
            
            beamAxis = (obj.nodes{end}.location - obj.nodes{1}.location)';
            beamAxis = beamAxis./norm(beamAxis);
            if dot(beamAxis, axis) ~= 0 
                error('Bending axis must be orthogonal to the beam axis');
            end
            currentBendingAxis = obj.ports{1}.orientation;
            alpha = acos(dot(currentBendingAxis, axis));
            R = round(vrrotvec2mat([beamAxis', alpha]), 12);     
            if sum(abs(R*currentBendingAxis - axis)) > 1e-9
                R = R';
            end
            
            for p=1:obj.n_ports
                obj.ports{p}.orientation = R*obj.ports{p}.orientation;
                obj.ports{p}.orientation(abs(obj.ports{p}.orientation) < 1e-9) = 0;
                obj.ports{p}.orientation = obj.ports{p}.orientation/norm(obj.ports{p}.orientation);
                
                if obj.ports{p}.orientation(2) == 1 && strcmp(obj.ports{p}.type, 'torque')
                    obj.ports{p}.orientation = -obj.ports{p}.orientation;
                end
            end
            
        end
    end
end