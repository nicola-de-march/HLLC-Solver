function FTV = TVSplitting(QL,QR)

global gamma g
global pdetype

switch pdetype
    case 0
        % =======================
        % SHALLOW WATERS
        % =======================
        % Taken from: A flux-vector splitting scheme for the shallow water 
        % equations extended to high-order on unstructured meshes, Toro et
        % al.
        hL = QL(1);            hR = QR(1);
        uL = QL(2)/hL;         uR = QR(2)/hR;
        psiL = QL(3)/hL;       psiR = QR(3)/hR;
        
        qL = hL*uL;
        qR = hR*uR;
        
        Qstar = 0.5 * (qL + qR) + (1/3)*sqrt(g)*(hL^(3/2) - hR^(3/2));
        hStar = (0.5*(hL^(3/2) + hR^(3/2)) - (3/(4*sqrt(g))) * (qR - qL));
        hStar = hStar^(2/3);

        if(Qstar>=0)
            F_a = [0; 
                   Qstar*uL;
                   Qstar*psiL];

        else
            F_a = [0; 
                   Qstar*uR;
                   Qstar*psiR];
        end
        F_p = [Qstar;
               0.5 * g* hStar^2;
               0];

        FTV = F_a + F_p;  % Combine the fluxes for the shallow water case
    case 1
        % =======================
        % EULER
        % =======================
        VL = PDECons2Prim(QL);
        VR = PDECons2Prim(QR);
        rhoL = VL(1);
        rhoR = VR(1);
        uL = VL(2);
        uR = VR(2);
        pL = VL(3);
        pR = VR(3);
        
        % PRESSURE SYSTEM:   
        % the solution for the pressure system (which is subsonic) is obtained via
        % linearised Riemann invariants.
        
        cL = sqrt(gamma*pL/rhoL);                                                   % speed of sound
        cR = sqrt(gamma*pR/rhoR);
        
        AL = sqrt(uL^2 + 4*cL^2);                                                   % for ideal gases
        AR = sqrt(uR^2 + 4*cR^2);
        CL = rhoL * (uL - AL);
        CR = rhoR * (uR + AR);
        
        ustar = (CR*uR - CL*uL)/(CR - CL) - 2*(pR - pL)/(CR - CL);                  % ustar_{i + 1/2}
        pstar = (CR*pL - CL*pR)/(CR - CL) + 0.5*(CR*CL/(CR - CL))*(uR - uL);        % pstar_{i + 1/2}
        
        if (ustar>=0)
        
            rhoK = rhoL;
        
        else
        
            rhoK = rhoR;
        
        end
        
        e = pstar / ( rhoK * (gamma-1) );
        Pflux = [0; pstar; ustar*(rhoK*e+pstar)];                                   % compute the pressure flux 
        
        
        %% ADVECTION SYSTEM: 
        % compute the advection flux based on the sign of ustar.
        
        if (ustar>=0)
            
            Aflux = ustar * [VL(1); VL(1)*VL(2); 0.5*VL(1)*VL(2)^2];
        
        else 
            
            Aflux = ustar * [VR(1); VR(1)*VR(2); 0.5*VR(1)*VR(2)^2];
        
        end
        
        
        %% Total TV numerical flux 
        
        FTV = Aflux + Pflux;


end