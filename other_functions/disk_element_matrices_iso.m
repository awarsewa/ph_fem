function [M, K] = disk_element_matrices_iso(nodes_xy, rho, E, nu, h)
% Check that we are given 4 corner nodes with x and y coordinats
if size(nodes_xy, 1) ~= 4 && size(nodes_xy, 2) ~= 2
    error('nodes_xy must be a 4x2 array of planar node coordinates');
end
% Check the orientation of the polygon
O = [1 nodes_xy(1, :); 1 nodes_xy(2, :); 1 nodes_xy(3,:)];
if sign(det(O)) <= 0
    error('please give the nodes in counter-clockwise order');
end

nodes_r = linspace(-1, 1, 2);
nodes_s = linspace(-1, 1, 2);

[~, phi_r] = lagrange_poly(nodes_r);
[~, dr_phi_r] = lagrange_poly_deriv(nodes_r, 1);

[~, phi_s] = lagrange_poly(nodes_s);
[~, ds_phi_s] = lagrange_poly_deriv(nodes_s, 1);


seq_r = [1 2 2 1];
seq_s = [1 1 2 2];
N = 4;
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

% Jacobian
Ji = cell(2);
Ji{1,1} = @(r,s) 0;
Ji{1,2} = @(r,s) 0;
Ji{2,1} = @(r,s) 0;
Ji{2,2} = @(r,s) 0;
for i = 1:4
    Ji{1,1} = @(r,s) Ji{1,1}(r, s) + dr_phi_rs{i}(r,s).*nodes_xy(i, 1); 
    Ji{1,2} = @(r,s) Ji{1,2}(r, s) + dr_phi_rs{i}(r,s).*nodes_xy(i, 2);
    Ji{2,1} = @(r,s) Ji{2,1}(r, s) + ds_phi_rs{i}(r,s).*nodes_xy(i, 1);
    Ji{2,2} = @(r,s) Ji{2,2}(r, s) + ds_phi_rs{i}(r,s).*nodes_xy(i, 2);
end

N_w = 2* (N - 2) - 1;
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
M = zeros(2*N);
for i=1:n_w
    for j =1:n_w
        Jk = zeros(2);
        Jk(1, 1) = J_i(i, j, 1);
        Jk(1, 2) = J_i(i, j, 2);
        Jk(2, 1) = J_i(i, j, 3);
        Jk(2, 2) = J_i(i, j, 4);
        N_ij = [phi_i(i,j,1), 0, phi_i(i,j,2), 0, phi_i(i,j,3), 0, phi_i(i, j, 4) 0;
                0, phi_i(i,j,1), 0, phi_i(i,j,2), 0, phi_i(i,j,3), 0, phi_i(i, j, 4)];

        M = M +  N_ij'*N_ij.*(det(Jk)*w(j)*w(i));
    end
end

M = rho.*h.*M;

D = [E*h/(1-nu^2), E*h*nu/(1-nu^2), 0;
     E*h*nu/(1-nu^2), E*h/(1-nu^2), 0;
     0, 0, E*h/(2+2*nu)];
 
K = zeros(2*N);
for i=1:n_w
    for j =1:n_w
        Jk = zeros(2);
        Jk(1, 1) = J_i(i, j, 1);
        Jk(1, 2) = J_i(i, j, 2);
        Jk(2, 1) = J_i(i, j, 3);
        Jk(2, 2) = J_i(i, j, 4);
        invJ = inv(Jk);
        B_ij = zeros(3, 2*N);
        for k = 1:N
            B_ij(1, 2*k-1) = dr_phi_i(i,j,k)*invJ(1,1)+ds_phi_i(i,j,k)*invJ(1,2);
            B_ij(2, 2*k) = dr_phi_i(i,j,k)*invJ(2,1)+ds_phi_i(i,j,k)*invJ(2,2); 
            B_ij(3, 2*k-1) =  B_ij(2, 2*k); 
            B_ij(3, 2*k) = B_ij(1, 2*k-1);
        end

        K = K +  B_ij'*D*B_ij.*(det(Jk)*w(j)*w(i));
    end
end
 
end