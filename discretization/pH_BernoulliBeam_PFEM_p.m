function [n, J, Q, G] = pH_BernoulliBeam_PFEM_p(N_p, N_q, myu, E, I, L) 
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
    
    % Cumpute int_Z phi_pzz * phi_q' dz --> D
    % First we need the Gauss-Legendre nodes + weights for N_p + N_q
    N_d = N_p + N_q - 4;
    n_d = ceil((N_d +1)/2); 
    [d, w_d] = legpts(n_d, interval);
    % Now compute the integral...
    D = zeros(N_p, N_q);
    phi_pd = deriv(phi_p, d, 2);
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
    
    %{
    % Input/output mappings (clamped-free boundary conditions)
    Tq1 = [     deriv(phi_p,interval(1),1);    % Ma
                deriv(phi_p,interval(2),1)];   % Mb
    Tq2 = [     phi_p(interval(1));            % Fa
                phi_p(interval(2))];           % Fb
                
    Bp = [ Tq1', Tq2' ];
    %}
    
    Bp = [  deriv(phi_p,linspace(0, L, max(N_q-3, 2))',1)', phi_p(linspace(0, L, max(N_q-3,2))')'];
    
    % Energy matrix
    Qq = E*I*M_q^-1;
    Qp = 1/myu*M_p^-1;
    Q = [Qp zeros(N_p, N_q);
        zeros(N_q, N_p) Qq]; 
    
    % Structure matrix
    J = [zeros(N_p, N_q) -D;
         D' zeros(N_q, N_p)];
    % Input\matrix + feedthrough matrix
    G = [Bp; zeros(N_p, size(Bp, 2))];
    
    decimalPrecision = 14;  
    J = round(J, decimalPrecision);
    Q = round(Q, decimalPrecision);
    G = round(G, decimalPrecision);
      
    % System order
	n = N_p+N_q;
end