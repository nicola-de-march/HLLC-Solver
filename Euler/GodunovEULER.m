
global g 
global pdetype rm_solver

% loop over all spatial control volumes 
Qnew = zeros(3,IMAX);
QGod = zeros(3, 1);
VL = zeros(3,IMAX);
for i=1:IMAX
    switch rm_solver
        case 0
            % Exact Riemann Solver
            s = 0.0;

            if (i==1)
                % Dirichlet BC on the left 
                % Godunov flux
                %
                VL = PDECons2Prim(Q(:, i));
                VR = PDECons2Prim(Q(:, i+1));

                [QGod(1), QGod(2), QGod(3)] = ExactRiemannEuler(VL(1), VR(1), ...
                                         VL(2), VR(2), ...
                                         VL(3), VR(3), ...
                                         s, gamma, gamma);
                QGod = PDEPrim2Cons(QGod);
                fp = f(QGod);
                [QGod(1), QGod(2), QGod(3)] = ExactRiemannEuler(V_init_L(1), VR(1), ...
                                         V_init_L(2), VR(2), ...
                                        V_init_L(3), VR(3), ...
                                         s, gamma, gamma);
                QGod = PDEPrim2Cons(QGod);
                fm = f(QGod);
            elseif (i==IMAX)
                 % Dirichlet BC on the right 
                % Godunov flux
                VL = PDECons2Prim(Q(:, i-1));
                VR = PDECons2Prim(Q(:, i));

                [QGod(1), QGod(2), QGod(3)] = ExactRiemannEuler(VL(1), V_init_R(1), ...
                                         VL(2), V_init_R(2), ...
                                         VL(3), V_init_R(3), ...
                                         s, gamma, gamma);

                QGod = PDEPrim2Cons(QGod);
                fp = f(QGod);
                [QGod(1), QGod(2), QGod(3)] = ExactRiemannEuler(VL(1), VR(1), ...
                                         VL(2), VR(2), ...
                                         VL(3), VR(3), ...
                                         s, gamma, gamma);
                QGod = PDEPrim2Cons(QGod);
                fm = f(QGod);
            else
                % Godunov flux
                VL = PDECons2Prim(Q(:, i));
                VR = PDECons2Prim(Q(:, i+1));
                [QGod(1), QGod(2), QGod(3)] = ExactRiemannEuler(VL(1), VR(1), ...
                                         VL(2), VR(2), ...
                                         VL(3), VR(3), ...
                                         s, gamma, gamma);
                QGod = PDEPrim2Cons(QGod);
                fp = f(QGod);

                VL = PDECons2Prim(Q(:, i-1));
                VR = PDECons2Prim(Q(:, i));
                [QGod(1), QGod(2), QGod(3)] = ExactRiemannEuler(VL(1), VR(1), ...
                                         VL(2), VR(2), ...
                                         VL(3), VR(3), ...
                                         s, gamma, gamma);
                QGod = PDEPrim2Cons(QGod);
                fm = f(QGod);
            end
        case 1
            % HLLC 
            if (i==1)
                % Dirichlet BC on the left 
                fp = HLLCSolver(Q(:,i),Q(:,i+1));
                fm = HLLCSolver(QL,Q(:,i));
            elseif (i==IMAX)
                 % Dirichlet BC on the right 
                fp = HLLCSolver(Q(:,i),QR);
                fm = HLLCSolver(Q(:,i-1),Q(:,i));
            else
                fp = HLLCSolver(Q(:,i),Q(:,i+1));
                fm = HLLCSolver(Q(:,i-1),Q(:,i));
            end
        case 2
            % Rusanov
            if (i==1)
                % Dirichlet BC on the left 
                fp = Rusanov(Q(:,i), Q(:,i+1));
                fm = Rusanov(QL, Q(:,i));
            elseif (i==IMAX)
                 % Dirichlet BC on the right 
                fp = Rusanov(Q(:,i), QR);
                fm = Rusanov(Q(:,i-1), Q(:,i));
            else
                fp = Rusanov(Q(:,i), Q(:,i+1));
                fm = Rusanov(Q(:,i-1), Q(:,i));
            end
        case 3
            % TV Splitting
            if (i==1)
                % Dirichlet BC on the left 
                fp = TVSplitting(Q(:,i), Q(:,i+1));
                fm = TVSplitting(QL, Q(:,i));
            elseif (i==IMAX)
                 % Dirichlet BC on the right 
                fp = TVSplitting(Q(:,i), QR);
                fm = TVSplitting(Q(:,i-1), Q(:,i));
            else
                fp = TVSplitting(Q(:,i), Q(:,i+1));
                fm = TVSplitting(Q(:,i-1), Q(:,i));
            end
    end
    % FV method
    Qnew(:,i) = Q(:,i) - dt/dx*(fp-fm);
end  
