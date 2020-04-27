function [n, J, Q, G] = pH_Link_PFEM_p(N_p, N_q, myu, E, A, L) 
    % Choose number of nodes for flow space
    interval = [0, L];
    nodes_p = linspace(0, L, N_p);
    nodes_q = linspace(0, L, N_q);
    
    % Construct polynomial bases (Lagrange polynomials)
    phi_p = lagrange_poly(nodes_p);
    phi_q = lagrange_poly(nodes_q);
    
    % For the weak form of the Stokes-Dirac structure choose the flow bases
    % as test functions v(z)

    % Cumpute int_Z phi_pz * phi_q' dz --> D
    % First we need the Gauss-Legendre nodes + weights for N_p + N_q
    N_d = N_p + N_q - 3;
    n_d = ceil((N_d +1)/2); 
    [d, w_d] = lgwt(n_d, interval(1), interval(2));
    % Now compute the integral...
    D = zeros(N_p, N_q);
    phi_pd = lagrange_poly_deriv(nodes_p, 1);
    phi_pd = phi_pd(d);
    phi_qd = phi_q(d);
    for i=1:N_p
        for j=1:N_q
            D(i,j) = w_d' * (phi_pd(:, i) .* phi_qd(:, j));
        end
    end 
    
    % Compute M_p
    N_mp = 2* N_p -2;
    n_mp = ceil((N_mp +1)/2); 
    [mp, w_mp] = lgwt(n_mp, interval(1), interval(2));
    M_p = zeros(N_p,N_p);
    phi_pm = phi_p(mp);
    for i=1:N_p
        for j=1:N_p
            M_p(i,j) = w_mp' * (phi_pm(:, i) .* phi_pm(:, j));
        end
    end 
    
    % Compute M_q
    N_mq = 2* N_q -2;
    n_mq = ceil((N_mq +1)/2); 
    [mq, w_mq] = lgwt(n_mq,  interval(1), interval(2));
    M_q = zeros(N_q,N_q);
    phi_qm = phi_p(mq);
    for i=1:N_q
        for j=1:N_q
            M_q(i,j) = w_mq' * (phi_qm(:, i) .* phi_qm(:, j));
        end
    end 
    Bp = phi_p(linspace(0, L, max(N_q,2))')';
    
    % Energy matrix
    Qq = E*A*M_q^-1;
    Qp = 1/myu*M_p^-1;
    Q = [Qp zeros(N_p, N_q);
        zeros(N_q, N_p) Qq]; 

    % Structure matrix
    J = [zeros(N_p, N_q) -D;
         D' zeros(N_q, N_p)];
    % Input matrix + feedthrough matrix
    G = [Bp; zeros(N_q, size(Bp, 2))];
    
    decimalPrecision = 14;
    
    % Substitute user-supplied values in calculated matrices
    J = round(J, decimalPrecision);
    Q = round(Q, decimalPrecision);
    G = round(G, decimalPrecision);
    
    n = N_p + N_q;
end