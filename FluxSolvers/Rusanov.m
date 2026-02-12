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
        pL = QL(3); % left pressure
        pR = QR(3); % right pressure
        
        fL = [QL(1)*uL; 
              QL(1)*uL^2 + pL; 
              QL(1)*uL*psiL];
        fR = [QR(1)*uR; 
              QR(1)*uR^2 + pR; 
              QR(1)*uR*psiR];
        
        % Calculate maximum wave speed for Euler equations
        smax = max(abs(uL) + sqrt(gamma * pL / QL(1)), abs(uR) + sqrt(gamma * pR / QR(1)));
        
        % Rusanov flux for Euler equations
        flux = 0.5 * (fL + fR) - 0.5 * smax * (QR - QL);
end
