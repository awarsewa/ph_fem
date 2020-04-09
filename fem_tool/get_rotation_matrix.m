function C = get_rotation_matrix(dim, elem_type, elem_nodes)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% Load type defines:
load elem_types type_*

% Calculate element vector:
elem_vec = elem_nodes(2, :) - elem_nodes(1, :);
L = sqrt(sum(elem_vec.*elem_vec));

switch dim
    % 2-dimensional case:
    case 2
        % angle of rotation:
        cos_alpha = sum(elem_vec.*[1,0]) / sqrt(sum(elem_vec.*elem_vec));
        if cos_alpha ~= 0
            sin_alpha = sign(cos_alpha) * sign(elem_vec(1)*elem_vec(2)) * sqrt(1 - cos_alpha^2);
        else
            sin_alpha = 1;
        end
        
        % Calculate rotation matrix for nodes:
        switch elem_type
            case type_link
                C = [cos_alpha, -sin_alpha; sin_alpha, cos_alpha];
            case type_beam
                C = [cos_alpha, -sin_alpha, 0; sin_alpha, cos_alpha, 0; 0, 0, 1];
        end
    
    % 3-dimensional case:
    case 3
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
        
        % Calculate rotation matrix for nodes:
        switch elem_type
            case type_link
                C = [n(1)*n(1)*(1-cos_alpha) +      cos_alpha, n(1)*n(2)*(1-cos_alpha) - n(3)*sin_alpha, n(1)*n(3)*(1-cos_alpha) + n(2)*sin_alpha; ...
                     n(2)*n(1)*(1-cos_alpha) + n(3)*sin_alpha, n(2)*n(2)*(1-cos_alpha) +      cos_alpha, n(2)*n(3)*(1-cos_alpha) - n(1)*sin_alpha; ...
                     n(3)*n(1)*(1-cos_alpha) - n(2)*sin_alpha, n(3)*n(2)*(1-cos_alpha) + n(1)*sin_alpha, n(3)*n(3)*(1-cos_alpha) +      cos_alpha];
            case type_beam
                C = [n(1)*n(1)*(1-cos_alpha) +      cos_alpha, n(1)*n(2)*(1-cos_alpha) - n(3)*sin_alpha, n(1)*n(3)*(1-cos_alpha) + n(2)*sin_alpha; ...
                     n(2)*n(1)*(1-cos_alpha) + n(3)*sin_alpha, n(2)*n(2)*(1-cos_alpha) +      cos_alpha, n(2)*n(3)*(1-cos_alpha) - n(1)*sin_alpha; ...
                     n(3)*n(1)*(1-cos_alpha) - n(2)*sin_alpha, n(3)*n(2)*(1-cos_alpha) + n(1)*sin_alpha, n(3)*n(3)*(1-cos_alpha) +      cos_alpha];
                C = [C, 0.*C; 0.*C, C];
        end        
        
end

% Build rotation matrix for element:
C = [C, 0.*C; 0.*C, C];

end

