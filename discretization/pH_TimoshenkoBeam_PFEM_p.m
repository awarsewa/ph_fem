function [n, J, Q, G] = pH_TimoshenkoBeam_PFEM_p(N_p, N_q, rho, E, G, A, I, kappa, L) 
    addpath '../.'
    % Interpolation interval
    interval = [0, L];

    % NOTE: For the folliwing, the Chebfun package is required!
    %       See http://www.chebfun.org/ for details.
    %       After downloading Chebfun, make sure to add it to the Matlab PATH

    % collocation points are the Gauss-Legendre nodes
    %z_p = legpts(N_p, interval);
    %z_q = legpts(N_q, interval);
    % Construct polynomial bases (Lagrange polynomial)
    %phi_p = chebfun.lagrange(z_p, interval);
    %phi_q = chebfun.lagrange(z_q, interval); 
    phi_p = chebfun.lagrange(linspace(0, L, N_p), interval);
    phi_q = chebfun.lagrange(linspace(0, L, N_q), interval);
    
    
    % Cumpute int_Z phi_pz * phi_q' dz --> D
    % First we need the Gauss-Legendre nodes + weights for N_p + N_q
    N_d = N_p + N_q - 3;
    n_d = ceil((N_d +1)/2); 
    [d, w_d] = legpts(n_d, interval);
    % Now compute the integral...
    D = zeros(N_p, N_q);
    phi_pd = deriv(phi_p, d, 1);
    phi_qd = phi_q(d);
    for i=1:N_p
        for j=1:N_q
            D(i,j) = w_d * (phi_pd(:, i) .* phi_qd(:, j));
        end
    end 
    
    % Compute M_p
    N_mp = 2* N_p -2;
    n_mp = ceil((N_mp +1)/2); 
    [mp, w_mp] = legpts(n_mp, interval);
    M_p = zeros(N_p,N_p);
    phi_pm = phi_p(mp);
    for i=1:N_p
        for j=1:N_p
            M_p(i,j) = w_mp * (phi_pm(:, i) .* phi_pm(:, j));
        end
    end 
    
    % Compute M_q
    N_mq = 2* N_q -2;
    n_mq = ceil((N_mq +1)/2); 
    [mq, w_mq] = legpts(n_mq, interval);
    M_q = zeros(N_q,N_q);
    phi_qm = phi_q(mq);
    for i=1:N_q
        for j=1:N_q
            M_q(i,j) = w_mq * (phi_qm(:, i) .* phi_qm(:, j));
        end
    end 
    
    % Input/output mappings (clamped-free boundary conditions)
    Tq = [     phi_p(interval(1));            % Ma/Fa
               phi_p(interval(2))];           % Mb/Fb
            
    Bp = [ Tq', zeros(N_p, 2);      % M
           zeros(N_p, 2), Tq'];     % F

    % Energy matrix
    Qq1 = E*I*M_q^-1;
    Qq2 = kappa*A*G*M_q^-1;
    Qp1 = 1/(rho*I)*M_p^-1;
    Qp2 = 1/(rho*A)*M_p^-1;
    Q = [Qp1 zeros(N_p, 2*N_q+N_p);
        zeros(N_p), Qp2, zeros(N_p, 2*N_q);
        zeros(N_q, 2*N_p) Qq1, zeros(N_q);
        zeros(N_q, 2*N_p+N_q), Qq2]; 

    % Structure matrix
    J = [zeros(N_p, 2*N_q) -D M_p;
         zeros(N_p, 2*N_q + N_p) -D;
         D' zeros(N_q, 2*N_p + N_q);
         -M_q D' zeros(N_q, 2*N_p)];
    % Input matrix + feedthrough matrix
    G = [Bp; zeros(2*N_q, 4)];
    
    decimalPrecision = 14;  
    J = round(J, decimalPrecision);
    Q = round(Q, decimalPrecision);
    G = round(G, decimalPrecision);
      
    % System order
	n = 2*(N_p+N_q);
end