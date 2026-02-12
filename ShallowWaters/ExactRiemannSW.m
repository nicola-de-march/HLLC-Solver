% Exact Riemann solver for shallow water equations (St. Venant)
% Input:
%   - QL, QR : left and right state vectors
%   - xi     : similarity parameter (xi=x/t)
% Output:
%   - Q      : solution Q=Q(xi)
function Q = ExactRiemannSW(QL,QR,xi)
    global g            % gravity constant
    hL = QL(1);         % left depth
    hR = QR(1);         % right depth
    uL = QL(2)/QL(1);   % left velocity
    uR = QR(2)/QR(1);   % right velocity
    psiL = QL(3)/QL(1); % left concentration
    psiR = QR(3)/QR(1); % right concentration
    hs = Newton(hL,hR,uL,uR); % compute h* using Newton
    us = uL - phi(hs,hL);     % compute u* by definition

    % "sampling" of the solution
    if(xi<=us)
        % left of us
        psi = psiL; 
        if(hs>hL)
            % left shock
            s = uL - sqrt(0.5*g*hs/hL*(hL+hs)); 
            if(xi<=s)
                h=hL;
                u=uL; 
            else
                h=hs;
                u=us; 
            end
        else
            % left rarefaction
            head = uL - sqrt(g*hL); 
            tail = us - sqrt(g*hs); 
            if(xi<=head) 
                % left
                h=hL; 
                u=uL; 
            elseif(xi>=tail)
                % right
                h=hs;
                u=us; 
            else
                % inside rarefaction
                h = ( (uL+2*sqrt(g*hL)-xi)/3 )^2 / g; 
                u = xi + sqrt(g*h); 
            end
        end
    else
        % on the right of us
        psi=psiR;     
        if(hs>hR)
            % right shock
            s = uR + sqrt(0.5*g*hs/hR*(hs+hR)); 
            if(xi<=s)
                h = hs; 
                u = us; 
            else
                h = hR;
                u = uR; 
            end
        else
            % right rarefaction
            tail = us + sqrt(g*hs); 
            head = uR + sqrt(g*hR); 
            if(xi<=tail)
                h=hs;
                u=us;
            elseif(xi>=head)
                h=hR;
                u=uR;
            else
                h = ( (xi-uR+2*sqrt(g*hR))/3 )^2 / g;
                u = xi - sqrt(g*h); 
            end
        end
    end
    % conservative variables vector
    Q = [h; h*u; h*psi]; 
end