% Function of the physical flux for SW
function pflux = f(Q)
    global g gamma
    global pdetype

    switch pdetype
        case 0
            % =========================
            % SHALLOW WATER 1D
            % =========================
            pflux = zeros(3,1);  % force Matlab to employ a column vector
            h = Q(1);
            if (h <= 0)
                h=0;
                u=0;
                psi=0;
            else
                u   = Q(2)/Q(1); 
                psi = Q(3)/Q(1); 
            end
            pflux(1) = h*u;                  % mass flux
            pflux(2) = h*u^2 + 0.5*g*h^2;    % momentum flux
            pflux(3) = h*u*psi;              % mass flux for the transported substance
        case 1
            % =========================
            % EULER 1D
            % =========================
            pflux = zeros(3,1);

            rho = Q(1);
            mom = Q(2);
            E   = Q(3);

            rho_eps = 1e-14;

            if rho < rho_eps
                u = 0;
                p = 0;
            else
                u = mom/rho;
                p = (gamma-1)*(E - 0.5*rho*u^2);
            end

            pflux(1) = mom;
            pflux(2) = mom*u + p;
            pflux(3) = u*(E + p);
    end
end