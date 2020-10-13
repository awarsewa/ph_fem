classdef PH_FluidicBeam < PH_NonlinearSystem
    properties(SetAccess = protected, GetAccess = public)
        N
        density
        youngsModulus
        shearModulus
        areaMomentOfInertia
        crossSectionalArea
        timoshenkoShearCoefficient
        length
        fluidBulkModulus
        pressureChambers
        alpha1
        alpha2
    end
    
    methods(Access = public)
        %              PH_NonlinearSystem(n, J, H, G, R)
        function obj = PH_FluidicBeam(N, rho, E, G, A, I, kappa, L, alpha1, alpha2, beta, pc_struct)
            [n_, J_, Q_, G_, phi_pd] = pH_TimoshenkoBeam_PFEM_p(N, N, rho, E, G, A, I, kappa, L);
            % Damping
            Q_1 = Q_(1:2*N, 1:2*N);
            Q_2 = Q_(2*N+1:end, 2*N+1:end);
            J_12 = J_(1:2*N, 2*N+1:end);
            J_21 = J_(2*N+1:end, 1:2*N);
            
            D = alpha1*(Q_1)^-1 - alpha2 * J_12*Q_2*J_21;
            D = D - 0.5*(D - D');
            R_ = casadi.MX.zeros(4*N, 4*N);
            R_(:,:) = blkdiag(D, zeros(2*N));
            
            n_pc = length(pc_struct);
            n_w = ceil((N-1)/2);
            g_pc = zeros(N, n_pc);
            x_ = casadi.MX.sym('x', n_ + 2*n_pc, 1);
            H_ = x_(1:n_)'*Q_*x_(1:n_);
            for k = 1:n_pc
                r = 0.5*(pc_struct(k).d^2 - pc_struct(k).c^2);
                [p, w] = lgwt(n_w, pc_struct(k).a, pc_struct(k).b);
                g_pc(:, k) = r.*phi_pd(p)'*w;
                
                J_ = [J_, [-g_pc(:, k), (beta/x_(n_+(k-1)*2+1))*g_pc(:, k); casadi.MX.zeros(n_-N+2*(k-1), 2)]; ...
                      [g_pc(:, k)'; -(beta/x_(n_+(k-1)*2+1))*g_pc(:, k)'], zeros(2, n_-N+2*k)];
                H_ = H_ + pc_struct(k).A*x_(n_+(k-1)*2+1)*(beta*(exp(x_(n_+2*k)/beta)-1)-x_(n_+2*k));
                G_ = [G_, zeros(n_+2*(k-1), 2); [zeros(2, 2*N+2*(k-1)), ...
                                                [0; beta/(x_(n_+2*(k-1)+1)*pc_struct(k).A)*pc_struct(k).k_v*sqrt(pc_struct(k).p_s - x_(n_+2*k))], ...
                                                [0; beta/(x_(n_+2*(k-1)+1)*pc_struct(k).A)*pc_struct(k).k_v*sqrt(x_(n_+2*k) - pc_struct(k).p_t)]]];
                R_ = blkdiag(R_, zeros(2));
            end
            n_ = n_ + 2*n_pc;
            obj = obj@PH_NonlinearSystem(n_, x_, J_, H_, G_, R_);
            
            obj.name = 'Timoshenko beam with integrated fluidic actuators';
            
            % Physical parameters
            obj.N = N;
            obj.youngsModulus = E;
            obj.shearModulus = G;
            obj.areaMomentOfInertia = I;
            obj.crossSectionalArea = A;
            obj.length = L;
            obj.density = rho;
            obj.timoshenkoShearCoefficient = kappa;
            obj.fluidBulkModulus = beta;
            obj.pressureChambers = pc_struct; 
            obj.alpha1 = alpha1;
            obj.alpha2 = alpha2;
            
            obj.n_nodes = N;
            node_pos = linspace(0, L, obj.n_nodes);
            for n = 1:obj.n_nodes
                obj.nodes{n} = PH_MechanicalNode([n, obj.n_nodes+n], 1, [node_pos(n) 0 0]);
            end
            
            obj.elements{1}.delete();
            %PH_Element(type, nodes, ports, varargin)
            obj.elements{1} = PH_Element('PH_FluidicBeam', 1:obj.n_nodes, 1:2*N+n_pc); 
            
            obj.n_ports = 2*N + 2*n_pc;
            for n = 1:obj.n_nodes
                %PH_MechanicalPort_boundary(type, nodes, orientation, element, IOPair, inputName, outputName) 
                obj.ports{n} = PH_MechanicalPort_boundary('torque', n, [0; 1; 0], 1, n, ['M' num2str(n) '_y'], ['w' num2str(n) '_y']);
                obj.ports{obj.n_nodes+n} = PH_MechanicalPort_boundary('force', n, [0; 0; 1], 1, obj.n_nodes+n, ['F' num2str(n) '_z'], ['v' num2str(n) '_z']);
            end
            for p = 1:n_pc
                obj.ports{2*N+2*(p-1)+1} = PH_HydraulicPort_external('valve', 2*N+2*(p-1)+1, ['x_v' num2str(p) ' >= 0'], ['h_' num2str(p)]);
                obj.ports{2*N+2*p} = PH_HydraulicPort_external('valve', 2*N+2*p, ['x_v' num2str(p) ' <= 0'], ['h_' num2str(p)]);
            end
        end
        
        function setPosition(obj, nodePositions)
            if ~ismatrix(nodePositions) || ...
                size(nodePositions, 1) ~= 2 || ...
                size(nodePositions, 2) ~= 3 
                error('NodePositions must be a 2x3 matrix');
            end
            dir = obj.nodes{2}.location - obj.nodes{1}.location;
            dir = dir./norm(dir);
            nDir = nodePositions(2,:) - nodePositions(1,:);
            nDir = nDir./norm(nDir);
            
            Tr = round(vrrotvec2mat(vrrotvec(dir,nDir)), 12);
            
            for p=1:2*obj.N
                obj.ports{p}.orientation = Tr*obj.ports{p}.orientation;
                obj.ports{p}.orientation(abs(obj.ports{p}.orientation) < 1e-9) = 0;
                obj.ports{p}.orientation = obj.ports{p}.orientation/norm(obj.ports{p}.orientation);
            end    

            for n=1:obj.n_nodes
                obj.nodes{n}.location = nodePositions(n,:);
            end
        end
    end
    
    methods(Access = protected)   
        function copyData(obj, cp)
            copyData@PH_NonlinearSystem(obj, cp);
        end
        function cp = copyElement(obj)
               % PH_FluidicBeam(N, rho, E, G, A, I, kappa, L, beta, pc_struct)
            cp = PH_FluidicBeam(obj.N, obj.density, obj.youngsModulus, obj.shearModulus, obj.crossSectionalArea, ... 
                                obj.secondMomentOfArea, obj.timoshenkoShearCoefficient, obj.length, ...
                                obj.alpha1, obj.alph2, obj.fluidBulkModulus, obj.pressureChambers);
            obj.copyData(cp);
        end
    end
end