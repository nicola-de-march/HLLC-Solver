% Function to compute the eigenvalues
function L = Lambda(Q)

    global g gamma
    global pdetype

    Q = Q(:);   % forza colonna

    switch pdetype

        % =========================
        % SHALLOW WATER 1D
        % =========================
        case 0

            h = Q(1);
            h_eps = 1e-14;
            if h < h_eps
                u = 0;
                c = 0;
            else
                u = Q(2)/h;
                c = sqrt(g*h);
            end

            L = [u-c;
                 u;
                 u+c];

        % =========================
        % EULER 1D
        % =========================
        case 1

            rho = Q(1);
            mom = Q(2);
            E   = Q(3);

            rho_eps = 1e-14;

            if rho < rho_eps
                u = 0;
                c = 0;
            else
                u = mom / rho;
                p = (gamma-1) * (E - 0.5*rho*u^2);
                c = sqrt(gamma * p / rho);
            end

            L = [u-c;
                 u;
                 u+c];
    end

end
