% Derivative of the function to be solved
function y=dfunc(hs,hL,hR,uL,uR) 
    eps = 1e-7; 
    % central finite difference approximation
    y = (func(hs+eps,hL,hR,uL,uR)-func(hs-eps,hL,hR,uL,uR))/(2*eps); 
end