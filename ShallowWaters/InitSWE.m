%% Define the Riemann problem
hL   = 0.15;                % left water depth
hR   = 0.05;                % right water depth
uL   = 0;                   % left velocity
uR   = 0;                   % right velocity
psiL = 1;                   % left scalar 
psiR = 0;                   % right scalar 

QL = [hL; hL*uL; hL*psiL];  % left state vector
QR = [hR; hR*uR; hR*psiR];  % right state vector

% init Q
Q = zeros(3,IMAX);
for i=1:IMAX
    if (x(i)<=0) 
        Q(:,i) = QL;
    else
        Q(:,i) = QR;
    end
end