function V = PDECons2Prim(Q)
% Compute primitivies from conservative variales
global gamma 
% Euler
V = Q;

V(1, :) = Q(1, :);
V(2, :) = Q(2, :) ./ Q(1, :);
V(3, :) = (gamma - 1) * (Q(3, :) - 0.5*Q(2, :).^2./Q(1, :));
