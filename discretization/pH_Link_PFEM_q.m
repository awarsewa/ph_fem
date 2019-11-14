function [n, J, Q, G] = pH_Link_PFEM_q(N_p, N_q, myu, E, A, L) 
    % Choose number of nodes for flow space
    interval = [0, L];
    
    % NOTE: For the following, the Chebfun package is required!
    %       See http://www.chebfun.org/ for details.
    %       After downloading Chebfun, make sure to add it to the Matlab PATH

    % collocation points are the Gauss-Legendre nodes
    z_p = legpts(N_p, interval);
    z_q = legpts(N_q, interval);
    % Construct polynomial bases (Lagrange polynomial)
    phi_p = chebfun.lagrange(z_p, interval);
    phi_q = chebfun.lagrange(z_q, interval); 

    % For the weak form of the Stokes-Dirac structure choose the flow bases
    % as test functions v(z)

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
    phi_qm = phi_p(mq);
    for i=1:N_q
        for j=1:N_q
            M_q(i,j) = w_mq * (phi_qm(:, i) .* phi_qm(:, j));
        end
    end 

    % Input/output mappings 
    % If we apply integration by parts to the impulse equations, our inputs 
    % are forces 
    Tq1 = phi_q(interval(1));   % v_a
    Tq2 = phi_q(interval(2));   % v_b       

    Bq = [ -Tq1', Tq2' ];

    % Energy matrix
    Qq = E*A*M_q^-1;
    Qp = 1/myu*M_p^-1;
    Q = [Qp zeros(N_p, N_q);
        zeros(N_q, N_p) Qq]; 

    % Structure matrix
    J = [zeros(N_p, N_q) D';
         -D zeros(N_q, N_p)];
    % Input matrix + feedthrough matrix
    G = [zeros(N_q, 2); Bq];
    
    decimalPrecision = 14;
    
    % Substitute user-supplied values in calculated matrices
    J = round(J, decimalPrecision);
    Q = round(Q, decimalPrecision);
    G = round(G, decimalPrecision);
    
    n = N_p + N_q;
end