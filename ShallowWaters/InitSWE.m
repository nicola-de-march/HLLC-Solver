%% Define the Riemann problem
hL   = 0.51;                % left water depth
hR   = 0.48;                % right water depth
uL   = 2.5;                   % left velocity
uR   = -5.8;                   % right velocity
psiL = 1;                   % left scalar 
psiR = 0;                   % right scalar 
xc = 0;

QL = [hL; hL*uL; hL*psiL];  % left state vector
QR = [hR; hR*uR; hR*psiR];  % right state vector

% init Q
Q = zeros(3,IMAX);
for i=1:IMAX
    if (x(i)<=xc) 
        Q(:,i) = QL;
    else
        Q(:,i) = QR;
    end
end