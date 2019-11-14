classdef PH_MechanicalPort < PH_Port 
    % Mechanical ports have an additional orientation in x/y/z coordinates
    properties (Access = public) 
        orientation % Orientation of this port in x/y/z coordinates
    end
    methods (Access = public)
        function obj = PH_MechanicalPort(scope, type, IOPair, inputType, outputType, node, orientation)
            if ~(strcmp(type, 'force/velocity') || strcmp(type, 'torque/angular velocity') ...
                    || strcmp(type, 'integrator'))
                error('type must be either ''force/velocity'' or ''torque/angular velocity'' or ''integrator''.');
            end
            
            if strcmp(inputType, outputType)
                error('inputType and outputType cannot be identical'); 
            end

            if strcmp(type, 'force/velocity')
                if ~(strcmp(inputType, 'force') || strcmp(inputType, 'velocity'))
                    error('inputType must be either ''force'' or ''velocity''');
                end
                if ~(strcmp(outputType, 'force') || strcmp(outputType, 'velocity'))
                    error('outputType must be either ''force'' or ''velocity''');
                end
            end
            
            if strcmp(type, 'torque/angular velocity')
                if ~(strcmp(inputType, 'torque') || strcmp(inputType, 'angular velocity'))
                    error('inputType must be either ''torque'' or ''angular velocity''');
                end
                if ~(strcmp(outputType, 'torque') || strcmp(outputType, 'angular velocity'))
                    error('outputType must be either ''torque'' or ''angular velocity''');
                end
            end
            
            obj = obj@PH_Port('mechanical', scope, type, IOPair, inputType, outputType, node);
            
            if nargin > 6
                if isempty(orientation) || ~iscolumn(orientation) || length(orientation)  ~= 3
                    error('orientation must be supplied as a 3x1 vector')
                end
                obj.orientation = orientation;
            else 
               obj.orientation = zeros(3, 1);                
            end
        end
        
        function obj = PH_Port(obj)
            obj = PH_Port('mechanical', obj.scope, obj.type, obj.IOPair, obj.inputType, obj.outputType, obj.node);
        end
    end
end