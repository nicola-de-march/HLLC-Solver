% Newton method, needed to compute h* 
function hs=Newton(hL,hR,uL,uR)
    global g 
    % initial guess
    %hs = 0.5*(hL+hR);  % arithmetic average (simpler)
    hs = ( sqrt(hL)+sqrt(hR) - (uR-uL)/2/sqrt(g) )^2/4; 
    tol = 1e-12;       % tolerance
    MaxNewton = 100;   % maximum number of iteractions
    for iNewton=1:MaxNewton
        gk = func(hs,hL,hR,uL,uR);  % function to solve
        res = abs(gk); 
        %disp(sprintf(' iNewton = %d, res = %e ', iNewton, res)) 
        if(res<tol)        
            break                   % we have found the solution so we exit
        end
        dg = dfunc(hs,hL,hR,uL,uR); % derivative of the funtion
        dh = -gk/dg;                % Newton step
        % line search globalization 
        delta = 1;                  % scale factor 0<delta<=1 to reduce Newton step, if necesary
        for inner=1:20         
            if( abs(func(hs+dh*delta,hL,hR,uL,uR)) >= res )
                % if the residual of the next iteration increases,
                % then the Newton step is reduced by 2 (factor:2)
                delta = 0.5*delta; 
            else 
                % residual does not increase => exit inner loop
                break
            end
        end
        % Update the solution. For delta=1, globalization method reduces to
        % standard Newton. 
        hs = hs + delta*dh; 
    end
end