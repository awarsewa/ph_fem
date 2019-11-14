function [n, J, Q, G, P, S, M] = pH_Link_PSM(N, mu_u, E_u, A_u, L_u) 
    % Choose number of nodes for flow space
    N_f = N;
    interval = [0, L_u];
    % Numner of nodes in the effort space
    N_e = N_f + 1;

    
    % NOTE: For the following, the Chebfun package is required!
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

    % Now compute the integral...
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

    % Same for .. int_Z psi * psi_z' dz
    N_d1 = (N_f-1) + (N_e - 2);
    n_d1 = ceil((N_d1+1)/2);
    [d1, w_d1] = legpts(n_d1, interval);
    phi_e_d1 = phi(d1);
    psi_ez = deriv(psi, d1, 1);


    D1_bar = zeros(N_f, N_e);
    for i=1:N_f
        for j=1:N_e
            D1_bar(i,j) = w_d1 * (phi_e_d1(:,i) .* psi_ez(:, j));
        end
    end

    syms myu E A
    assume([myu E A], 'real');
    assume([myu E A], 'positive');

    % Input/output mappings (clamped-free boundary conditions)
    Tp1 = psi(interval(1));     % v_a
    Tq1 = zeros(1, N_e);        
    Tq2 = psi(interval(2));     % F_b
    Tp2 = zeros(1, N_e);

    Sq1 = -psi(interval(1));    % -F_a
    Sp1 = zeros(1, N_e);
    Sp2 = psi(interval(2));     % v_b
    Sq2 = zeros(1, N_e); 

    T = [ Tq1 Tp1;
          Tq2 Tp2];
    S = [ Sq1 Sp1;
          Sq2 Sp2];

    D1 = M_phi^(-1) * D1_bar;

    E_ = [  zeros(N_f, N_e)     D1;
            D1                  zeros(N_f,N_e);
            S];

    F_ = [  M_psiphi'            zeros(N_f, N_e);
            zeros(N_f, N_e)      M_psiphi';
            T];   

    sys = E_ * F_^-1;

    % Energy matrix
    Q1 = E*A*M_phi;
    Q2 = 1/myu*M_phi;
    Q = [Q1 zeros(N_f);
        zeros(N_f) Q2]; 

    % Structure matrix
    J = sys(1:2*N_f, 1:2*N_f);

    % Input matrix + feedthrough matrix
    G = sys(1:2*N_f, 2*N_f+1:end);
    D = sys(2*N_f+1:end, 2*N_f+1:end);

    % Submatrices needed for damping assignment
    %J12 = J(1:N_f,N_f+1:end);
    %J21 = J(N_f+1:end, 1:N_f);

    G_t = sys(2*N_f+1:end, 1:2*N_f);

    M =  D;
    S = sym(zeros(size(M)));
    S(2,2) = M(2, 2);
    M(2,2) = 0;

    % Do not add dissipation yet
    %R = [   zeros(N_f,2*N_f);
             %zeros(N_f), b0*Q2^(-1)-b1*J21*Q1*J12];
    
    %{
    A_ = J - R;
    % Make J skew-symmetric
    R = sym(zeros(2*N_f,2*N_f)); 
    for i=1:N_f
        for j=1:N_f
            if double(subs(A_(N_f+i,N_f+j), [E A myu], [E_u A_u mu_u])) < double(subs(A_(N_f+j,N_f+i), [E A myu], [E_u A_u mu_u]))
                R(N_f+i,N_f+j) = (-A_(N_f+i, N_f+j) - A_(N_f+j, N_f+i))/2;
                R(N_f+j,N_f+i) = (-A_(N_f+i, N_f+j) - A_(N_f+j, N_f+i))/2;
                A_(N_f+i,N_f+j) = A_(N_f+i,N_f+j) + R(N_f+i,N_f+j);
                A_(N_f+j,N_f+i) = A_(N_f+j,N_f+i) + R(N_f+j,N_f+i);
            else
                if(i == j)
                    R(N_f+i,N_f+j) = -A_(N_f+i, N_f+j); 
                    A_(N_f+i,N_f+j) = 0; 
                end
            end
        end
    end
    J = A_;    
    %}
    
    
    % Factor a matrix P out of B and G such that G_star equals B_star'
    if any(any(G' ~= G_t))
        P = sym(zeros(2,2*N_f));
        G = G';
        for i=1:2
            for j=1:2*N_f
                P(i,j) = (G_t(i,j)-G(i,j))/2;
            end
        end
        G = (G + P)';
        P = P';
    else
        P = zeros(2*N_f, 2);
    end
    
    decimalPrecision = 12; 
    
    % Substitute user-supplied values in calculated matrices
    J = round(double(subs(J, [myu, E, A], [mu_u, E_u, A_u])),decimalPrecision);
    Q = round(double(subs(Q, [myu, E, A], [mu_u, E_u, A_u])),decimalPrecision);
    G = round(double(subs(G, [myu, E, A], [mu_u, E_u, A_u])),decimalPrecision);
    P = round(double(subs(P, [myu, E, A], [mu_u, E_u, A_u])),decimalPrecision);
    S = round(double(subs(S, [myu, E, A], [mu_u, E_u, A_u])),decimalPrecision);
    M = round(double(subs(M, [myu, E, A], [mu_u, E_u, A_u])),decimalPrecision);
    
    n = 2*N_f;
end