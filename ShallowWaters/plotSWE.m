function plotSWE(Q, x, t)

% Plot the numerical solution (left part in blue, right part in yellow, as in the laboratory experiment)
psi = Q(3,:)./Q(1,:);
u   = Q(2,:)./Q(1,:);
leftc  = psi>0.5;
rightc = psi<=0.5; 
fig1=figure(1);
fig1.Position=[500   540   560   420]; 
subplot(3,1,1)
hold off 
plot(x(leftc),Q(1,leftc),'bo')                                                                  % h on the left of the contact 
hold on 
plot(x(rightc),Q(1,rightc),Color=[0.9290 0.6940 0.1250],Marker='o',LineStyle='none')            % h on the right of the contact 
title(sprintf('Time = %f',t))
xlabel('x')
ylabel('h')
subplot(3,1,2)
hold off 
plot(x(leftc),u(leftc),'bo')                                                                    % u on the left of the contact 
hold on 
plot(x(rightc),u(rightc),Color=[0.9290 0.6940 0.1250],Marker='o',LineStyle='none')              % u on the right of the contact 
xlabel('x')
ylabel('u')
subplot(3,1,3)
hold off 
plot(x(leftc),psi(leftc),'bo')                                                                  % psi on the left of the contact 
hold on 
plot(x(rightc),psi(rightc),Color=[0.9290 0.6940 0.1250],Marker='o',LineStyle='none')            % psi on the right of the contact 
xlabel('x')
ylabel('psi')