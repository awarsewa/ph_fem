function [n, J, Q, G, P, S, M] = pH_BernoulliBeam_PSM(N, myu, E, I, L) 
    addpath '../.'
    % Choose number of nodes for flow space
    N_f = N;
    % Interpolation interval
    interval = [0, L];
    
    % Numner of nodes in the effort space
    N_e = N_f + 2;

    % NOTE: For the folliwing, the Chebfun package is required!
    %       See http://www.chebfun.org/ for details.
    %       After downloading Chebfun, make sure to add it to the Matlab PATH

    % collocation points are the Gauss-Legendre nodes
    z_f = legpts(N_f, interval);
    z_e = legpts(N_e, interval);

    % Construct polynomial bases (Lagrange polynomial)
    phi = chebfun.lagrange(z_f, interval);
    psi = chebfun.lagrange(z_e, interval);

    % For the weak form of the Stokes-Dirac structure choose the effort bases
    % as test functions v(z)

    % Cumpute int_Z psi * phi' dz
    % First we need the Gauss-Legendre nodes + weights for N_f + N_e
    N_m = N_f + N_e -2;
    n_m = ceil((N_m +1)/2); 
    [m, w_m] = legpts(n_m, interval);

    % Now compute the integral...phi
    M_psiphi = zeros(N_e, N_f);
    phi_m = phi(m);
    psi_m = psi(m);

    for i=1:N_e
        for j=1:N_f
            M_psiphi(i,j) = w_m * (psi_m(:, i) .* phi_m(:, j));
        end
    end 

    N_phi = 2* N_f -2;
    n_phi = ceil((N_phi +1)/2); 
    [n, w_n] = legpts(n_phi, interval);
    
    M_phi = zeros(N_f,N_f);
    phi_n = phi(n);
    for i=1:N_f
        for j=1:N_f
            M_phi(i,j) = w_n * (phi_n(:, i) .* phi_n(:, j));
        end
    end 

    % Same for .. int_Z psi * psi_zz' dz
    N_d2 = (N_f-1) + (N_e - 3);
    n_d2 = ceil((N_d2+1)/2);
    [d2, w_d2] = legpts(n_d2, interval);
    phi_e_d2 = phi(d2);
    psi_ezz = deriv(psi, d2, 2);

    D2_bar = zeros(N_f, N_e);
    for i=1:N_f
        for j=1:N_e
            D2_bar(i,j) = w_d2 * (phi_e_d2(:,i) .* psi_ezz(:, j));
        end
    end

    % Input/output mappings (clamped-free boundary conditions)
    
    Tp = [     psi(interval(1));               % w_a
               -deriv(psi,interval(1),1)];      % v_a
    Tq = [     -deriv(psi,interval(2),1);       % F_b
               psi(interval(2))];              % M_b
    Sq = [     deriv(psi,interval(1),1);       % F_a
               psi(interval(1))];              % M_a
    Sp = [     psi(interval(2))                % v_b
               deriv(psi,interval(2),1)];      % w_b
            
    T_ = [ Tp zeros(2,N_e);
          zeros(2,N_e) Tq];
    S_ = [ zeros(2,N_e) Sq;
          Sp zeros(2, N_e)];
          
    D2 = M_phi^(-1) * D2_bar;

    E_ = [  zeros(N_f, N_e)     -D2;
            D2                 zeros(N_f,N_e);
            S_];

    F_ = [   M_psiphi'           zeros(N_f, N_e);
             zeros(N_f, N_e)      M_psiphi';
             T_];   

    sys = E_ * F_^-1;

    % Energy matrix
    Qq = E*I*M_phi;
    Qp = 1/myu * M_phi;
    Q = [Qp zeros(N_f);
        zeros(N_f) Qq]; 

    % Structure matrix
    J = sys(1:2*N_f, 1:2*N_f);

    % Input matrix + feedthrough matrix
    G = sys(1:2*N_f, 2*N_f+1:end);
    D = sys(2*N_f+1:end, 2*N_f+1:end);
    
    G_t = sys(2*N_f+1:end, 1:2*N_f);

    M =  D;
    S = zeros(size(M));
    S(3:4,3:4) = M(3:4, 3:4);
    M(3:4,3:4) = 0;

    R = zeros(2*N_f);
    
    A = J - R;
    % Make J skew-symmetric
    R = zeros(2*N_f,2*N_f); 
    for i=1:N_f
        for j=1:N_f
            if A(N_f+i,N_f+j) < A(N_f+j,N_f+i)
                R(N_f+i,N_f+j) = (-A(N_f+i, N_f+j) - A(N_f+j, N_f+i))/2;
                R(N_f+j,N_f+i) = (-A(N_f+i, N_f+j) - A(N_f+j, N_f+i))/2;
                A(N_f+i,N_f+j) = A(N_f+i,N_f+j) + R(N_f+i,N_f+j);
                A(N_f+j,N_f+i) = A(N_f+j,N_f+i) + R(N_f+j,N_f+i);
            else
                if(i == j)
                    R(N_f+i,N_f+j) = -A(N_f+i, N_f+j); 
                    A(N_f+i,N_f+j) = 0; 
                end
            end
        end
    end
    J = A;    

    % Factor a matrix P out of B and G such that G_star equals B_star'
    if any(any(G' ~= G_t))
        P = zeros(4,2*N_f);
        G = G';
        for i=1:4
            for j=1:2*N_f
                P(i,j) = (G_t(i,j)-G(i,j))/2;
            end
        end
        G = (G + P)';
        P = P';
    else
        P = zeros(2*N_f, 4);
    end

    decimalPrecision = 11; 

    % Substitute user-supplied values in calculated matrices
    J = round(J, decimalPrecision);
    Q = round(Q, decimalPrecision);
    G = round(G, decimalPrecision);

    P = round(P, decimalPrecision);
    S = round(S, decimalPrecision);
    M = round(M, decimalPrecision);
    
    n = 2*N_f;
end