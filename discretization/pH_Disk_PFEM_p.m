function [n, J, Q, G, nodes_xy] = pH_Disk_PFEM_p(order, nodes_xy, rho, E, nu, h)
% Check that we are given 4 corner nodes with x and y coordinats
if size(nodes_xy, 1) ~= 4 && size(nodes_xy, 2) ~= 2
    error('nodes_xy must be a 4x2 array of planar node coordinates');
end
% Check the orientation of the polygon
O = [1 nodes_xy(1, :); 1 nodes_xy(2, :); 1 nodes_xy(3,:)];
if sign(det(O)) <= 0
    error('please give the nodes in counter-clockwise order');
end
if order < 1
    error('order of polynomials must be at least 1');
end

% number of nodes 
N = (order+1)^2;
seq_r = [1 order+1 order+1 1];
seq_s = [1 1 order+1 order+1];
if order > 1
    for i=1:order-1
        nodes_xy(4*i+1, :) = nodes_xy(1, :) + i/order*(nodes_xy(2, :) - nodes_xy(1, :));
        nodes_xy(4*i+2, :) = nodes_xy(2, :) + i/order*(nodes_xy(3, :) - nodes_xy(2, :));
        nodes_xy(4*i+3, :) = nodes_xy(3, :) - i/order*(nodes_xy(3, :) - nodes_xy(4, :));
        nodes_xy(4*i+4, :) = nodes_xy(4, :) - i/order*(nodes_xy(4, :) - nodes_xy(1, :));
        seq_r(4*i+1:4*i+4) = [1+i order+1 order+1-i 1];
        seq_s(4*i+1:4*i+4) = [1 i+1 order+1 order+1-i];
    end
    % number of inner rectangles
    N_in = floor((order+1)/2-1);
    for i = 1:N_in
        inner(1, :) = nodes_xy(1, :) + i/order*(nodes_xy(2, :) - nodes_xy(1, :)) + i/order*(nodes_xy(4, :) - nodes_xy(1, :));
        inner(2, :) = nodes_xy(2, :) - i/order*(nodes_xy(2, :) - nodes_xy(1, :)) + i/order*(nodes_xy(3, :) - nodes_xy(2, :));
        inner(3, :) = nodes_xy(3, :) - i/order*(nodes_xy(3, :) - nodes_xy(4, :)) - i/order*(nodes_xy(3, :) - nodes_xy(2, :));
        inner(4, :) = nodes_xy(4, :) - i/order*(nodes_xy(4, :) - nodes_xy(1, :)) + i/order*(nodes_xy(3, :) - nodes_xy(4, :));
        seq_r_inner = [1+i order+1-i order+1-i 1+i];
        seq_s_inner = [1+i 1+i order+1-i order+1-i];
        for j=1:order-2-i
            inner(4*j+1, :) = inner(1, :) + j/(order-2*i)*(inner(2, :) - inner(1, :));
            inner(4*j+2, :) = inner(2, :) + j/(order-2*i)*(inner(3, :) - inner(2, :));
            inner(4*j+3, :) = inner(3, :) - j/(order-2*i)*(inner(3, :) - inner(4, :));
            inner(4*j+4, :) = inner(4, :) - j/(order-2*i)*(inner(4, :) - inner(1, :));
            seq_r_inner(4*j+1:4*j+4) = [1+i+j order+1-i order+1-i-j 1+i];
            seq_s_inner(4*j+1:4*j+4) = [1+i i+1+j order+1-i order+1-i-j];
        end
        nodes_xy = [nodes_xy; inner]; 
        seq_r = [seq_r, seq_r_inner];
        seq_s = [seq_s, seq_s_inner];
    end
    % single node at the center
    if mod(order+1, 2)
        nodes_xy(end+1, :) =  nodes_xy(1,:) + 0.5*(nodes_xy(4,:) + 0.5*(nodes_xy(3,:) - nodes_xy(4,:)) - nodes_xy(1,:) + 0.5*(nodes_xy(2,:) - nodes_xy(1,:)));
        seq_r(end+1) = floor((order+1)/2)+1;
        seq_s(end+1) = floor((order+1)/2)+1;
    end
end

nodes_r = linspace(-1, 1, order+1);
nodes_s = linspace(-1, 1, order+1);

[~, phi_r] = lagrange_poly(nodes_r);
[~, dr_phi_r] = lagrange_poly_deriv(nodes_r, 1);

[~, phi_s] = lagrange_poly(nodes_s);
[~, ds_phi_s] = lagrange_poly_deriv(nodes_s, 1);

nodes_rs = zeros(N, 2); 
for i = 1:N
    nodes_rs(i, 1) = nodes_r(seq_r(i));
    nodes_rs(i, 2) = nodes_s(seq_s(i));
end

% build required functions of r and s and their derivatives
phi_rs = cell(1, N);
dr_phi_rs = cell(1, N);
ds_phi_rs = cell(1, N);
for i = 1:N
    phi_rs{i}       = @(r, s) phi_r{seq_r(i)}(r) .* phi_s{seq_s(i)}(s);
    dr_phi_rs{i}    = @(r, s) dr_phi_r{seq_r(i)}(r) .* phi_s{seq_s(i)}(s);
    ds_phi_rs{i}    = @(r, s) phi_r{seq_r(i)}(r) .* ds_phi_s{seq_s(i)}(s);
end

[~, Ni] = lagrange_poly([-1, 1]);
[~, dN] = lagrange_poly_deriv([-1, 1], 1);

Nr = cell(4,1);
Nr{1} = @(r, s) dN{1}(r).*Ni{1}(s);
Nr{2} = @(r, s) dN{2}(r).*Ni{1}(s);
Nr{3} = @(r, s) dN{2}(r).*Ni{2}(s);
Nr{4} = @(r, s) dN{1}(r).*Ni{2}(s);

Ns = cell(4,1);
Ns{1} = @(r, s) Ni{1}(r).*dN{1}(s);
Ns{2} = @(r, s) Ni{2}(r).*dN{1}(s);
Ns{3} = @(r, s) Ni{2}(r).*dN{2}(s);
Ns{4} = @(r, s) Ni{1}(r).*dN{2}(s);

% Jacobian
Ji = cell(2);
Ji{1,1} = @(r,s) 0;
Ji{1,2} = @(r,s) 0;
Ji{2,1} = @(r,s) 0;
Ji{2,2} = @(r,s) 0;
for i = 1:4
    Ji{1,1} = @(r,s) Ji{1,1}(r, s) + Nr{i}(r,s).*nodes_xy(i, 1); 
    Ji{1,2} = @(r,s) Ji{1,2}(r, s) + Nr{i}(r,s).*nodes_xy(i, 2);
    Ji{2,1} = @(r,s) Ji{2,1}(r, s) + Ns{i}(r,s).*nodes_xy(i, 1);
    Ji{2,2} = @(r,s) Ji{2,2}(r, s) + Ns{i}(r,s).*nodes_xy(i, 2);
end

N_w = 2*(order+1) - 1;
n_w = ceil((N_w +1)/2); 
[pw, w] = lgwt(n_w, -1, 1);

J_i = zeros(n_w, n_w, 4);
J_i(:,:,1) = Ji{1,1}(pw, pw');
J_i(:,:,2) = Ji{1,2}(pw, pw');
J_i(:,:,3) = Ji{2,1}(pw, pw');
J_i(:,:,4) = Ji{2,2}(pw, pw');

phi_i = zeros(n_w, n_w, N);
dr_phi_i = zeros(n_w, n_w, N);
ds_phi_i = zeros(n_w, n_w, N);
for i = 1:N
    phi_i(:,:,i) = phi_rs{i}(pw, pw');
    dr_phi_i(:,:,i) = dr_phi_rs{i}(pw, pw');
    ds_phi_i(:,:,i) = ds_phi_rs{i}(pw, pw');
end

% Gauss-Legendre Integration (siehe IBB FE-Skript)
M = zeros(N);
for i=1:N
    for j=1:N
        for k=1:n_w
            Jk = zeros(n_w, 4);
            Jk(:, 1) = J_i(k, :, 1);
            Jk(:, 2) = J_i(k, :, 2);
            Jk(:, 3) = J_i(k, :, 3);
            Jk(:, 4) = J_i(k, :, 4);
            j_m = zeros(n_w,1);
            for p=1:n_w
                j_m(p) = det([Jk(p, 1), Jk(p, 2); Jk(p, 3), Jk(p, 4)]);
            end
            
            f_m = j_m.*phi_i(k,:,i)'.*phi_i(k,:,j)';
            M(i, j) = M(i, j) + w'*f_m*w(k);
        end
    end
end

D_x = zeros(N);
for i=1:N
    for j=1:N
        for k=1:n_w
            Jk = zeros(n_w, 4);
            Jk(:, 1) = J_i(k, :, 1);
            Jk(:, 2) = J_i(k, :, 2);
            Jk(:, 3) = J_i(k, :, 3);
            Jk(:, 4) = J_i(k, :, 4);
            
            f_dx = zeros(n_w,1);
            for p = 1:n_w
                J = [Jk(p, 1), Jk(p, 2); ...
                     Jk(p, 3), Jk(p, 4)];
                invJ = inv(J);
                detJ = det(J);
                dx_phi = invJ(1,1)*dr_phi_i(k,p,i) + ...
                         invJ(1,2)*ds_phi_i(k,p,i);
                f_dx(p) = detJ.*dx_phi*phi_i(k,p,j);
            end
            
            D_x(i, j) = D_x(i,j) + w'*f_dx*w(k);
        end
    end
end

D_y = zeros(N);
for i=1:N
    for j=1:N
        for k=1:n_w
            Jk = zeros(n_w, 4);
            Jk(:, 1) = J_i(k, :, 1);
            Jk(:, 2) = J_i(k, :, 2);
            Jk(:, 3) = J_i(k, :, 3);
            Jk(:, 4) = J_i(k, :, 4);

            f_dy = zeros(n_w,1);
            for p=1:n_w
                J = [Jk(p, 1), Jk(p, 2); ...
                     Jk(p, 3), Jk(p, 4)];
                invJ = inv(J);
                detJ = det(J);
                ddy_phi = invJ(2,1)*dr_phi_i(k,p,i) + ...
                          invJ(2,2)*ds_phi_i(k,p,i);
                f_dy(p) = detJ.*ddy_phi*phi_i(k,p,j);
            end
            
            D_y(i, j) = D_y(i,j) + w'*f_dy*w(k);
        end
    end
end

J = [zeros(N, 2*N), -D_x, zeros(N), -D_y; ... 
     zeros(N, 3*N), -D_y, -D_x; ...
     D_x', zeros(N, 4*N); ... 
     zeros(N), D_y', zeros(N, 3*N); ...
     D_y', D_x', zeros(N, 3*N)];
 
Q = [1/(rho*h)*inv(M), zeros(N, 4*N); ...
     zeros(N), 1/(rho*h)*inv(M), zeros(N, 3*N); ...
     zeros(N, 2*N), E*h/(1-nu^2)*inv(M), E*h*nu/(1-nu^2)*inv(M), zeros(N);
     zeros(N, 2*N), E*h*nu/(1-nu^2)*inv(M), E*h/(1-nu^2)*inv(M), zeros(N);
     zeros(N, 4*N), E*h/(2+2*nu)*inv(M)];


% Boundary
G = [eye(N), zeros(N); ...
    zeros(N),  eye(N); ...
    zeros(3*N, 2*N)];

n = 5*N;

end