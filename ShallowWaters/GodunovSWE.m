
global g 
global pdetype rm_solver

% loop over all spatial control volumes 
Qnew = zeros(3,IMAX);
    
for i=1:IMAX
    switch rm_solver
        case 0
            % Exact Riemann Solver
            if (i==1)
                % Dirichlet BC on the left 
                % Godunov flux
                QGod = ExactRiemannSW(Q(:,i),Q(:,i+1),0);
                fp = f(QGod);
                QGod = ExactRiemannSW(QL,Q(:,i),0);
                fm = f(QGod);
            elseif (i==IMAX)
                 % Dirichlet BC on the right 
                % Godunov flux
                QGod = ExactRiemannSW(Q(:,i),QR,0);
                fp = f(QGod);
                QGod = ExactRiemannSW(Q(:,i-1),Q(:,i),0);
                fm = f(QGod);
            else
                % Godunov flux
                QGod = ExactRiemannSW(Q(:,i),Q(:,i+1),0);
                fp = f(QGod);
                QGod = ExactRiemannSW(Q(:,i-1),Q(:,i),0);
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
                fp = Rusanov(Q(:,i),Q(:,i+1));
                fm = Rusanov(QL,Q(:,i));
            elseif (i==IMAX)
                 % Dirichlet BC on the right 
                fp = Rusanov(Q(:,i),QR);
                fm = Rusanov(Q(:,i-1),Q(:,i));
            else
                fp = Rusanov(Q(:,i),Q(:,i+1));
                fm = Rusanov(Q(:,i-1),Q(:,i));
            end
        case 3
            if (i==1)
                % Dirichlet BC on the left 
                fp = TVSplitting(Q(:,i),Q(:,i+1));
                fm = TVSplitting(QL,Q(:,i));
            elseif (i==IMAX)
                 % Dirichlet BC on the right 
                fp = TVSplitting(Q(:,i),QR);
                fm = TVSplitting(Q(:,i-1),Q(:,i));
            else
                fp = TVSplitting(Q(:,i),Q(:,i+1));
                fm = TVSplitting(Q(:,i-1),Q(:,i));
            end
    end
    % FV method
    Qnew(:,i) = Q(:,i) - dt/dx*(fp-fm);
end  
