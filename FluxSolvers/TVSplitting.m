function FTV = TVSplitting(QL,QR)

global gamma
global pdetype

switch pdetype
    case 0
        % =======================
        % SHALLOW WATERS
        % =======================
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