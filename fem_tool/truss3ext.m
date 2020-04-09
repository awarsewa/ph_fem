function [M, K, B, C, dof_mapping] = truss3ext(nodes, elems, props, global_dof_order, output)
%TRUSS3 Assemble 3D-model of truss structure for given
%geometry and material properties
%   [M, K, B, C] = TRUSS3(nodes, elems, elem_types, props, global_dof_order, output)
%   outputs the mass, stiffness, input and output matrix for a 2D truss structure made up
%   of n nodes and e elements.
%   
%   nodes is an n x 4 vector, where the first column contains the node number,
%   the second to fourth column the 3D-coordinates (x,y,z) of the nodes.
%   
%   elems is an e x 6 vector, where the first column
%   contains the element number, the second and third column contain the node
%   numbers of the two nodes connected by the element, and the fourth column
%   contains the element property group to which the element belongs. The
%   fifth column is the element type number, with the following mapping
%     1 -- link
%     2 -- beam
%   The sixth column is the actuation type number, with the following mapping:
%     0 -- no actuation
%     1 -- parallel actuation
%     2 -- serial actuation
%   
%   elem_types is a t x 2 cell array, where the first colum is the element type
%   number, and the second is a string with the type name. Currently only
%   'beam' and 'link' are supported types.
%
%   props is a p x 5 vector, where p is the number of different property groups. 
%   The first column of props contains the property group number, refered
%   to in the fourth column of elems. The second to fifth column of props
%   contain the physical parameters: Young's modulus, density, cross-sectional area
%   and yield strength.
%
%   global_dof_order is a 1 by n vector containing the order in which the
%   local nodal degrees of freedom map to the global degrees of freedom.
%   
%   output defines the calculated C-matrix. The following outputs are possible:
%     'disp': yields the displacements of the nodal DOFs
%     'fabs': yields the element forces (compression negative)
%     'frel': yields the element forces relative to the respective yield strength
%     'leng': yields the change in element length
%

%% Element Matrices:
    % local element and local coordinates defined as:
    %   y
    %   ^
    %   |         P1             P2
    %   |         |--------------|
    % z o---> x   |--> x1        |--> x2
    %
    
    % Initialize all local DOFs to zero:
    loc_dofs = zeros(length(nodes(:,1)), 1);
    
    for ii = 1:size(elems,1)
        prop_index = find(props(:,1) == elems(ii,4));
        node_index = [find(nodes(:,1) == elems(ii,2)), find(nodes(:,1) == elems(ii,3))];
        
        % Element Vector ist Vector from P1 to P2:  P1 ---> P2
        elem_vec = nodes(node_index(2), 2:4) - nodes(node_index(1), 2:4);
        elem_props = props(prop_index, 2:end);
        prop.E = elem_props(1);
        prop.rho = elem_props(2);
        prop.A = elem_props(3);
        prop.Iy = elem_props(4);
        prop.Iz = elem_props(5);
        prop.G = elem_props(6);
        prop.It = elem_props(7);
        L = sqrt(sum(elem_vec.*elem_vec));

        % Rotation matrix parameters:
        % normalized vector of rotation:
        n = cross([1;0;0], elem_vec(:));
        if sum(n) ~= 0
            n = n / sqrt(sum(n.*n));
        end
        % angle of rotation:
        cos_alpha = sum(elem_vec.*[1,0,0]) / L;
        if cos_alpha == 0
            sin_alpha = 1;
        else
            sin_alpha = sin(acos(cos_alpha));
        end
        
        % Get Mass, Stiffness and Input matrix according to element type:
        [M_elem{ii}, K_elem{ii}, B_elem{ii}, L_elem{ii}] = get_element_matrix_ext(3, elems(ii, 5:6), ([prop.E, prop.rho, prop.A, prop.Iy, prop.Iz, prop.G, prop.It, L]));
        
        % Get transformation to global coordinates for local element DOFs:
        C_elem{ii} = get_rotation_matrix(3, elems(ii,5), [nodes(node_index(1), 2:4); nodes(node_index(2), 2:4)]);
        
        % element force output:
        F_elem{ii} = K_elem{ii}(1,:);
        
        % local degrees of freedom (link: 3, beam: 6):
        loc_dofs(elems(ii,2)) = max(loc_dofs(elems(ii,2)), length(L_elem{ii})/2);
        loc_dofs(elems(ii,3)) = max(loc_dofs(elems(ii,3)), length(L_elem{ii})/2);
        
    end
    %% Assemble:

    n = sum(loc_dofs(:));
    dof_mapping = zeros(n, 2);
    
    K = spalloc(n,n,ceil(0.1*n^2));
    M = K;
    B = spalloc(n, sum(elems(:,5) == 3), ceil(0.1*n*sum(elems(:,5) == 3)));
    
    C_fabs = sparse(zeros(size(elems,1), n));
    C_sabs = sparse(zeros(size(elems,1), n));
    C_srel = sparse(zeros(size(elems,1), n));
    C_leng = sparse(zeros(size(elems,1), n));

    nb_act = 0;
    
    for ii = 1:size(elems,1)
        % find starting indices for element nodes in global DOF-vector:
        node_index = [find(nodes(:,1) == elems(ii,3)), find(nodes(:,1) == elems(ii,2))];
        global_node_index = [find(global_dof_order == nodes(node_index(1),1)), find(global_dof_order == nodes(node_index(2),1))];
        dof_start = [sum(loc_dofs(1:(global_node_index(1)-1)))+1, sum(loc_dofs(1:(global_node_index(2)-1)))+1];
        
        switch elems(ii,5)
            case 1
                tmp = [dof_start(1):(dof_start(1)+2), dof_start(2):(dof_start(2)+2)];
                dof_mapping(tmp,:) = [[ones(3,1)*elems(ii,3); ones(3,1)*elems(ii,2)], [1:3, 1:3]'];
            case 2
                tmp = [dof_start(1):(dof_start(1)+5), dof_start(2):(dof_start(2)+5)];
                dof_mapping(tmp,:) = [[ones(6,1)*elems(ii,3); ones(6,1)*elems(ii,2)], [1:6, 1:6]'];
            otherwise
                disp('Only "link", "beam" or "actuator" elements are supported !');
                return;
        end
        
        %C_loc = [C_elem{ii}, zeros(size(C_elem{ii})); zeros(size(C_elem{ii})), C_elem{ii}];
        
        K_loc = C_elem{ii}*K_elem{ii}*C_elem{ii}';
        M_loc = C_elem{ii}*M_elem{ii}*C_elem{ii}';
        F_loc = F_elem{ii}*C_elem{ii}';

        K(tmp, tmp) = K(tmp, tmp) + K_loc;
        M(tmp, tmp) = M(tmp, tmp) + M_loc;
        
        % Elementkraft:
        C_fabs(ii, tmp) = F_loc;
        % Elementspannung:
        C_sabs(ii, tmp) = C_fabs(ii, tmp)./props(elems(ii, 4), 4);
        % Elementausnutzung:
		C_srel(ii, tmp) = C_sabs(ii, tmp)./props(elems(ii, 4), 9);
        
        if elems(ii,6) == 1
            % Add actuator force:
            nb_act = nb_act + 1;
            B(tmp, nb_act) = C_elem{ii}*B_elem{ii};
            C_fabs(ii,:) = [];
            C_sabs(ii,:) = [];
            C_srel(ii,:) = [];
        end
        
        % Elementlänge:
		C_leng(ii, tmp) = L_elem{ii}*C_elem{ii}';
    end
	
    C_disp = eye(size(K));
    
    switch output
        case 'disp'
            C = C_disp;
        case 'fabs'
            C = C_fabs;
        case 'sabs'
            C = C_sabs;
        case 'srel'
            C = C_srel;
        case 'leng'
            C = C_leng;
        otherwise
            C = C_disp;
            disp(['Warning: Unrecognized parameter for ''output'': ', output, newline, 'Using ''disp'' instead!']);
    end
    
end