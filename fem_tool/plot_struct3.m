function plot_struct3(N, E, ET, H, draw_resp, limits)
%PLOT_STRUCT3 Plot the structure of a 3D truss building
%    PLOT_STRUCT3(N, E, ET) generates a 3D-plot showing a truss structure
%    described by the node list N, the element list E and the element types
%    listed in ET.
%    For n number of nodes, N is an n-by-4 matrix, with the first column
%    containing the node number, the second to fourth column containing the
%    [X,Y,Z]-coordinates of the node. For m number of elements, E is an
%    m-by-3 matrix containing the element number in the first column and
%    the number of start and end node in the second and third column for each
%    element, respectively.
%
%    PLOT_STRUCT3D(N, E, H) uses the figure handle H for plotting.
%    
%    PLOT_STRUCT3D(N, E, H, limits) uses the values in the row vector limits as
%    axis limits, where limits = [x_lim, y_lim] with x_lim and y_lim row vectors
%    containing the individual axes limits for the x- and y-axis, respectively.
%

    %%
    % number of elements:
    n = size(E, 1);

    % find limits:
    if ~exist('limits', 'var')
        x_lim = [min(N(:,2)), max(N(:,2))];
        x_lim = x_lim + [-1,1]*(x_lim(2) - x_lim(1))*0.8;
        y_lim = [min(N(:,3)), max(N(:,3))];
        y_lim = y_lim + [-1,1]*(y_lim(2) - y_lim(1))*0.1;
		z_lim = [min(N(:,4)), max(N(:,4))];
		z_lim = z_lim + [-1,1]*(z_lim(2) - z_lim(1))*0.1;
        limits = [x_lim, y_lim, z_lim];
    end
    % sort nodes ascending:
    % [~,n_ind] = sort(nodes(:,1));
    % nodes = nodes(n_ind,:);
    
    if ~exist('draw_resp', 'var')
        draw_resp = 1;
    end
    
    if ~exist('H', 'var')
        figure;
        H = axes;
    end

    longest_elem = 0;
    shortest_elem = Inf;
    for ii = 1:n
        ka = E(ii, 2);    % starting node
        ke = E(ii, 3);    %   ending node
        ind_elem = [find(N(:,1)==ka,1), find(N(:,1)==ke,1)];
        longest_elem = max(longest_elem, norm(N(ind_elem(1),2:3) - N(ind_elem(2),2:3)));
        shortest_elem = min(shortest_elem, norm(N(ind_elem(1),2:3) - N(ind_elem(2),2:3)));
    end
    
    for ii = 1:n
        ka = E(ii, 2);    % starting node
        ke = E(ii, 3);    %   ending node
        ind_elem = [find(N(:,1)==ka,1), find(N(:,1)==ke,1)];
        % if required nodes are in node list continue with plotting in black:    
        if length(ind_elem) == 2
            switch E(ii,5)
%                 case 'actuator'
%                     plot3(H, N(ind_elem,2), N(ind_elem,3), N(ind_elem,4), 'k', 'LineWidth', 1, 'LineStyle', '--');
                case 2
                    plot3(H, N(ind_elem,2), N(ind_elem,3), N(ind_elem,4), 'k', 'LineWidth', 2, 'LineStyle', '-');
                case 1
                    plot3(H, N(ind_elem,2), N(ind_elem,3), N(ind_elem,4), 'k', 'LineWidth', 1, 'LineStyle', '-');
            end
%             plot3(H, N(ind_elem,2), N(ind_elem,3), N(ind_elem,4), 'k');
            hold on;
            %plot3(H, N(ind_elem,2), N(ind_elem,3), N(ind_elem,4), 'ko', 'MarkerSize', 8);
            % Plotting with opacity according to element length:
            % plot(H, N(ind_elem,2), N(ind_elem,3), 'Color', [1, 0, 0, ((norm(N(ind_elem(1),2:3) - N(ind_elem(2),2:3))-shortest_elem)/((longest_elem == shortest_elem) + (longest_elem - shortest_elem)))]*0.7 + 0.3);
            if (mod(ii, 100) == 0) && (draw_resp ~= 0)
                drawnow;
            end
        else
            error('Node number from E not found in N. Check consistency of E and N.');
        end
    end

    %plot(x_lim, [0, 0], 'k');
    hold off;
    axis equal;
    xlim(limits(1:2));
    ylim(limits(3:4));
	zlim(limits(5:6));
end