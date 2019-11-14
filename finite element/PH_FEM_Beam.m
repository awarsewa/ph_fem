classdef PH_FEM_Beam < PH_FEM_element
    properties(Access = protected)
        N_b
        N_l
        N_t
        massPerUnitLength
        youngsModulus
        shearModulus
        areaMomentOfInertia_y
        areaMomentOfInertia_z
        polarMomentOfInertia
        crossSectionalArea
        torsionConstant
        length
    end
    methods(Access = public)
        function obj = PH_FEM_Beam(N_b, N_l, N_t, myu, E, A, G, I_x, I_y, I_z, L, kappa)
                % A complete beam consists of 2 Euler-Bernoulli beams, a
                % link for normal strain and a Saint Venant torsion element
                beam_type = 'Euler-Bernoulli';
                % If kappa is specified, the beam is of Timoshenko type!
                if nargin > 11
                    beam_type = 'Timoshenko';
                end
                
                if strcmp(beam_type, 'Timoshenko')
                    % PH_FEM_Timoshenko(N, rho, E, G, A, I, kappa, L)
                    bending_y = PH_FEM_Timoshenko(N_b, myu/A, E, G, A, I_y, kappa, L);
                    bending_z = PH_FEM_Timoshenko(N_b, myu/A, E, G, A, I_z, kappa, L);
                else
                    % PH_FEM_Bernoulli(N, myu, E, I, L, boundary)
                    bending_y = PH_FEM_Bernoulli(N_b, myu, E, I_y, L, 'free-free');
                    bending_z = PH_FEM_Bernoulli(N_b, myu, E, I_z, L, 'free-free');
                end
                
                % PH_FEM_Link(N, myu, E, A, L, boundary)
                normal  = PH_FEM_Link(N_l, myu, E, A, L, 'free-free');
                % PH_FEM_SaintVenant(N, I_p, G_t, I_t, L)
                torsion = PH_FEM_SaintVenant(N_t, I_y+I_z, G, I_x, L);
                
                bending_y.ports{1}.orientation = [0; 0; 1];
                bending_y.ports{2}.orientation = [0; 0; 1]; 
                bending_y.ports{3}.orientation = [0; 1; 0];
                bending_y.ports{4}.orientation = [0; 1; 0];
                bending_y.inputNames{1} = 'Ma_z';
                bending_y.inputNames{2} = 'Mb_z';
                bending_y.inputNames{3} = 'Fa_y';
                bending_y.inputNames{4} = 'Fb_y';
                bending_y.outputNames{1} = 'wa_z';
                bending_y.outputNames{2} = 'wb_z';
                bending_y.outputNames{3} = 'va_y';
                bending_y.outputNames{4} = 'vb_y';
                    
                obj = obj@PH_FEM_element(0, [], zeros(0), zeros(0), zeros(0));
                obj.add(bending_y);
                obj.add(bending_z);
                obj.add(torsion);
                obj.add(normal);

                obj.name = 'Port-Hamiltonian beam element';
                
                % Physical parameters
                obj.N_b = N_b;
                obj.N_l = N_l;
                obj.N_t = N_t;
                obj.massPerUnitLength = myu;
                obj.youngsModulus = E;
                obj.shearModulus = G;
                obj.areaMomentOfInertia_y = I_y;
                obj.areaMomentOfInertia_z = I_z;
                obj.polarMomentOfInertia = I_y + I_z;
                obj.crossSectionalArea = A;
                obj.torsionConstant = I_x;
                obj.length = L;
                
                obj.elements{1}.type = 'PH_FEM_Beam';
        end
    end
    
    
    
    methods(Access = protected)
        function copyData(obj, cp)
            copyData@PH_LinearSystem(obj, cp);
        end
         
        function cp = copyElement(obj)
            % PH_FEM_Beam(N_b, N_l, N_t, myu, E, A, G, I_x, I_y, I_z, L)
            cp = PH_FEM_Beam(obj.N_b, obj.N_l, obj.N_t, obj.massPerUnitLength, obj.youngsModulus, ...
                             obj.crossSectionalArea, obj.shearModulus, obj.torsionConstant, ..., 
                             obj.areaMomentOfInertia_y, obj.areaMomentOfInertia_z, obj.length);
            obj.copyData(cp);
        end
    end
end