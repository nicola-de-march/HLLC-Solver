function F = HLLCSolver(qL, qR)
% qL,qR: 3x1 conserved [h; hu; hpsi]
% F:   3x1 flux at interface
global g gamma
global pdetype

switch pdetype
    case 0
        % Shallow waters
        hmin = 1e-8;
        hL = max(qL(1), hmin); hR = max(qR(1), hmin);
        uL = qL(2)/hL;         uR = qR(2)/hR;
        psiL = qL(3)/hL;       psiR = qR(3)/hR;
        
        pL = 0.5*g*hL^2; pR = 0.5*g*hR^2;
        
        % flux vectors (3 components)
        fL = [qL(2); qL(2)*uL + pL; qL(2)*psiL];
        fR = [qR(2); qR(2)*uR + pR; qR(2)*psiR];
        
        % wave speed estimates
        cL = sqrt(g*hL); cR = sqrt(g*hR);
        SL = min(uL - cL, uR - cR);
        SR = max(uL + cL, uR + cR);
        
        if SL >= SR
            % fallback HLL average (3 comp.)
            F = (SR*fL - SL*fR + SR*SL*(qR - qL)) / (SR - SL);
            return
        end
        
        % contact speed S*
        Sstar = (SR*qR(2) - SL*qL(2) + pL - pR) / (SR*hR - SL*hL);
        
        % star states for h and hu
        hLstar = hL * (SL - uL) / (SL - Sstar);
        hRstar = hR * (SR - uR) / (SR - Sstar);
        qLstar = [hLstar; hLstar * Sstar; hLstar * psiL]; % third comp built from psiL
        qRstar = [hRstar; hRstar * Sstar; hRstar * psiR];
        
        % HLLC flux selection (3 comp.)
        if SL >= 0
            F = fL;
        elseif SL <= 0 && Sstar >= 0
            F = fL + SL*(qLstar - qL);
        elseif Sstar <= 0 && SR >= 0
            F = fR + SR*(qRstar - qR);
        else
            F = fR;
        end
    case 1 % Euler
        % Compute pressure and flux for Euler equations
        rhoL = qL(1); rhoR = qR(1);
        uL = qL(2)/rhoL;
        uR = qR(2)/rhoR;
        EL = qL(3); ER = qR(3);
        
        pL = (gamma-1)*(EL - 0.5*rhoL*uL^2);
        pR = (gamma-1)*(ER - 0.5*rhoR*uR^2);
        
        pL = max(pL,1e-10);
        pR = max(pR,1e-10);
        
        fL = [rhoL*uL;
              rhoL*uL^2 + pL;
              uL*(EL+pL)];
        
        fR = [rhoR*uR;
              rhoR*uR^2 + pR;
              uR*(ER+pR)];
        
        cL = sqrt(gamma*pL/rhoL);
        cR = sqrt(gamma*pR/rhoR);
        
        SL = min(uL-cL,uR-cR);
        SR = max(uL+cL,uR+cR);
        
        Sstar = (pR-pL + rhoL*uL*(SL-uL) - rhoR*uR*(SR-uR)) / ...
                (rhoL*(SL-uL) - rhoR*(SR-uR));
        
        rhoLstar = rhoL*(SL-uL)/(SL-Sstar);
        rhoRstar = rhoR*(SR-uR)/(SR-Sstar);
        
        ELstar = ((SL-uL)*EL - pL*uL + pL*Sstar)/(SL-Sstar);
        ERstar = ((SR-uR)*ER - pR*uR + pR*Sstar)/(SR-Sstar);
        
        qLstar = [rhoLstar; rhoLstar*Sstar; ELstar];
        qRstar = [rhoRstar; rhoRstar*Sstar; ERstar];
        
        if SL >= 0
            F = fL;
        elseif Sstar >= 0
            F = fL + SL*(qLstar-qL);
        elseif SR > 0
            F = fR + SR*(qRstar-qR);
        else
            F = fR;
        end

end

