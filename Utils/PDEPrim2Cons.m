% Compute conservative variables from primitive 
function Q = PDEPrim2Cons(V)
    global pdetype gamma 

    switch(pdetype)
        case 0 % Euler
            Q(1) = V(1); % density
            Q(2) = V(1) * V(2); % Momentum
            Q(3) = V(3) / (gamma - 1.0) + 0.5 * V(1) * V(2)^2;
        case 1 % SW with transport
            Q(1) = V(1);
            Q(2) = V(1) * V(2.3);
        case 3
            Q(1) = V(1);
            Q(2) = V(1) * V(2);
            Q(3) = V(3);
    end
end