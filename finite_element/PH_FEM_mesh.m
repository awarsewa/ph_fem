classdef PH_FEM_mesh < PH_LinearSystem
    properties(SetAccess = protected, GetAccess = public) 
        nodeTable
        elementTable    % array of elements and their connection to nodes
        elementTypes    % cell array containing the different element types
        attribs         % cell array of element attributes
    end
    
    methods(Access = public)
        function obj = PH_FEM_mesh(nodeArray, elementArray, elements, attributes)
            % Error checking 
            n_n = size(nodeArray, 1);
            n_e = size(elementArray, 1);
                        
            if n_n < 1 || n_e < 1
                error('At least one node/element has to be given');
            end
            
            if size(nodeArray, 2) ~= 3
                error('Node array must be an nx3 matrix');
            end
            
            if size(elementArray, 2) < 2
                error('Element array must have at least 2 entries per element. One node and the element type');
            end
            
            if ~isa(elements, 'cell')
                error('elements must be given as cell array of elements');
            end
            
            if nargin > 3 && ~isempty(attributes)
                if ~isa(attributes, 'cell') || size(attributes, 1) ~= n_e
                    error('attributes must be given as a cell array');
                end
            else
                attributes = cell(n_e, 1);
            end
            
            obj = obj@PH_LinearSystem(0, zeros(0), zeros(0), zeros(0));
            
            
            for n=1:n_n
                obj.nodes{n} = PH_MechanicalNode([], [], nodeArray(n, :));
            end
            obj.n_nodes = n_n;
            
            
            for i=1:n_e
                current = copy(elements{elementArray(i, end)});
                %current.setPosition(nodeArray(elementArray(i, 1:current.n_nodes), :));
                nodes = elementArray(i,1:end-1);
                nodes(nodes == 0) = [];
                current.setPosition(nodeArray(nodes, :));
                current.setAttributes(attributes(i, :));
                
                for k=1:length(current.ports)
                    current.ports{k}.inputName = ['e' num2str(i) '.' current.ports{k}.inputName];
                    current.ports{k}.outputName = ['e' num2str(i) '.' current.ports{k}.outputName];
                end
                % Concatenate the current element to the mesh
                obj.add(current);
                current.delete();
            end
  
            obj.elements{1}.delete();
            obj.elements(1) = [];
            obj.n_elements = obj.n_elements - 1;
            for n = 1:obj.n_nodes
                if ~isempty(obj.nodes{n}.elements)
                    obj.nodes{n}.elements = obj.nodes{n}.elements - 1; 
                end
            end
            
            % Set class properties
            obj.nodeTable = nodeArray;
            obj.elementTable = elementArray;
            obj.elementTypes = elements;
            obj.attribs = attributes;
        end
        
        % Modified add function that removes duplicate nodes
        function add(obj, system)
            obj.add@PH_LinearSystem(system);
            if isa(system, 'PH_FEM_mesh')  
                sys_elems = system.elementTable;
                sys_elems(:, 1:end-1) = sys_elems(:, 1:end-1) + size(obj.nodeTable, 1);
                sys_elems(:, end) = sys_elems(:, end) + size(obj.elementTypes, 2);
                obj.nodeTable = [obj.nodeTable; system.nodeTable];
                obj.elementTable = [obj.elementTable; sys_elems];
                obj.elementTypes = [obj.elementTypes, system.elementTypes];
                obj.attribs = [obj.attribs; system.attribs];
            end
        end 
        
        function sys = getSubsystem(obj, x_p, x_q)
            sys = obj.getSubsystem@PH_LinearSystem(x_p, x_q);
            % Remove nodes/elements that do not belong to the subsystem
            for n = size(sys.nodeTable, 1):-1:1
                if ~sys.getNodeAtLocation(sys.nodeTable(n,:))
                    for e = size(sys.elementTable, 1):-1:1
                        if any(sys.elementTable(e, 1:end-1) == n)
                            sys.elementTable(e, :) = [];
                            sys.attribs(e) = [];
                        end
                    end
                    sys.nodeTable(n,:) = [];
                end
            end
        end
        
        function pHLin = PH_LinearSystem(obj)
            pHLin = PH_LinearSystem(obj.n, obj.J, obj.Q, obj.G, ...
                                 obj.R, obj.P, obj.S, obj.M, obj.B, obj.x_p, obj.x_q, obj.C);
            pHLin.C_u = obj.C_u;
            pHLin.C_y = obj.C_y;
            pHLin.C_e = obj.C_e;
            pHLin.C_f = obj.C_f;
            pHLin.n_c = size(pHLin.C_u, 1) + size(pHLin.B, 2);

            pHLin.n_elements = obj.n_elements;
            pHLin.elements = cell(obj.n_elements, 1);
            for e=1:obj.n_elements
                pHLin.elements{e} = copy(obj.elements{e});
            end

            pHLin.n_nodes = obj.n_nodes;
            pHLin.nodes = cell(obj.n_nodes, 1);
            for n=1:obj.n_nodes
                pHLin.nodes{n} = copy(obj.nodes{n});
            end

            pHLin.n_ports = obj.n_ports;
            pHLin.ports = cell(obj.n_ports, 1);
            for p=1:obj.n_ports
                pHLin.ports{p} = copy(obj.ports{p});
            end
        end 
    end
    
    methods(Access = protected)        
        function copyData(obj, cp)
            for e = 1:cp.n_elements
                cp.elements{e}.delete();
            end
            for p = 1:cp.n_ports
                cp.ports{p}.delete(); 
            end
            for n = 1:cp.n_nodes
                cp.nodes{n}.delete();
            end
            copyData@PH_LinearSystem(obj, cp);
            cp.J = obj.J;
            cp.R = obj.R;
            cp.Q = obj.Q; 
            cp.G = obj.G;
            cp.P = obj.P;
            cp.S = obj.S;
            cp.M = obj.M;
            cp.C = obj.C;
            cp.x_q = obj.x_q;
            cp.x_p = obj.x_p;
            cp.n = obj.n;
            cp.n_u = obj.n_u;
            cp.C = obj.C;
            cp.B = obj.B;
        end
        function cp = copyElement(obj)           
            cp = PH_FEM_mesh(obj.nodeTable, obj.elementTable, obj.elementTypes, obj.attribs);
            obj.copyData(cp);
        end
    end
end