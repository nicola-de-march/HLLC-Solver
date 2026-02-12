% function containing the Rankine-Hugoniot relation, 
% Riemann invariants and entropy condition di Lax 
function y=phi(hs,hLR)
    global g 
    % Lax entropy condition
    if(hs>hLR) 
        % shock (Rankine-Hugoniot conditions)
        y = sqrt(0.5*g*(hs+hLR) / (hs*hLR)) * (hs-hLR);  
    else
        % rarefaction (Riemann invariants) 
        y = 2*sqrt(g)*( sqrt(hs) - sqrt(hLR) ); 
    end
end