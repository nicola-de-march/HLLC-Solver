%% Define the Riemann problem
global gamma

rhoL = 1;                % left density
rhoR = 0.125;            % right density
uL   = 0;                % left velocity
uR   = 0;                % right velocity
pL   = 1;                % left pressure 
pR   = 0.1;              % right pressure 
xc = 0;

V_init_L = [1; 0; 1];
V_init_R = [0.125; 0; 0.1];

QL = [rhoL; 
      rhoL*uL; 
      pL / (gamma - 1.0) + 0.5 * rhoL * uL^2];  % left state vector
QR = [rhoR; 
      rhoR*uR; 
      pR / (gamma - 1.0) + 0.5 * rhoR * uR^2];  % left state vector

% init Q
Q = zeros(3,IMAX);
for i=1:IMAX
    if (x(i)<=xc) 
        Q(:,i) = QL;
    else
        Q(:,i) = QR;
    end
end