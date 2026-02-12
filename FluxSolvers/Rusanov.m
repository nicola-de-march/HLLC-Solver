function flux = Rusanov(QL, QR)
% Rusanov flux for shallow water (supports [h;hu] or [h;hu;hpsi])
global g gamma
global pdetype

switch pdetype
    case 0
        % =========================
        % SHALLOW WATER 1D
        % =========================

        QL = QL(:); QR = QR(:);
        
        % primitives with protection
        hL = QL(1);         % left depth
        hR = QR(1);         % right depth
        uL = QL(2)/QL(1);   % left velocity
        uR = QR(2)/QR(1);   % right velocity
        psiL = QL(3)/QL(1); % left concentration
        psiR = QR(3)/QR(1); % right concentration
        
        % flux vectors
        pL = 0.5 * g * hL^2;
        pR = 0.5 * g * hR^2;
        
        % ---- Wave speeds ----
        cL = sqrt(g * hL);
        cR = sqrt(g * hR);
        
        fL = [hL*uL; 
              hL*uL^2 + pL; 
              hL*uL*psiL];
        fR = [hR*uR; 
              hR*uR^2 + pR; 
              hR*uR*psiR];
        
        % scalar dissipation speed (max char. speed)
        smax = max( abs(uL) + cL, abs(uR) + cR );
        
        % Rusanov flux
        flux = 0.5*(fL + fR) - 0.5*smax*(QR - QL);
    case 1
        % =========================
        % EULER 1D
        % =========================
        % Define flux vectors for Euler equations
        % Primitive variables
        VL = PDECons2Prim(QL);
        VR = PDECons2Prim(QR);

        rhoL = VL(1); 
        uL   = VL(2);
        pL   = VL(3);
        rhoR = VR(1); 
        uR   = VR(2);
        pR   = VR(3);
        
        EL = QL(3);
        ER = QR(3);
        
        % Physical fluxes
        fL = [
            rhoL*uL;
            rhoL*uL^2 + pL;
            uL*(EL + pL)
        ];
        
        fR = [
            rhoR*uR;
            rhoR*uR^2 + pR;
            uR*(ER + pR)
        ];
        
        % Sound speeds
        cL = sqrt(gamma * pL / rhoL);
        cR = sqrt(gamma * pR / rhoR);
        
        % Max wave speed
        smax = max(abs(uL) + cL, abs(uR) + cR);
        
        % Rusanov flux
        flux = 0.5*(fL + fR) - 0.5*smax*(QR - QL);

end
