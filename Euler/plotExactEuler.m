fig1=figure(1);
% Compute solution with the exact Riemann solver on a finer mesh 
EMAX = 10*IMAX; 
xe = linspace(xL,xR,EMAX); 
Ve = zeros(3,EMAX);
for i= 1:EMAX
    xi = xe(i)/time;
    [Ve(1,i), Ve(2, i), Ve(3, i)] = ExactRiemannEuler(V_init_L(1), V_init_R(1), ...
                             V_init_L(2), V_init_R(2), ...
                             V_init_L(3), V_init_R(3), ...
                             (xe(i)-xc)/time, gamma, gamma);
end

rho   = Ve(1, :); 
u     = Ve(2, :); 
p     = Ve(3, :);

subplot(3,1,1)
hold on
plot(xe, rho, 'b-')  


subplot(3,1,2)
hold on
plot(xe, u, 'b-')    


subplot(3,1,3)
hold on
plot(xe, p, 'b-')   

drawnow