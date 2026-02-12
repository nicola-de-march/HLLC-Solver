function plotEuler(Q, x, t)

% Converte Q in variabili primitive
V = PDECons2Prim(Q);  % V = [rho; u; p]

rho = V(1,:);
u   = V(2,:);
p   = V(3,:);

fig1 = figure(1);
fig1.Position = [500 540 560 420];

subplot(3,1,1)
cla
plot(x, rho, Color=[0.9290 0 0.1250], Marker='o', MarkerSize=3, LineStyle='none')
title(sprintf('Time = %f', t))
xlabel('x')
ylabel('\rho')

subplot(3,1,2)
cla
plot(x, u, Color=[0.9290 0 0.1250], Marker='o', MarkerSize=3, LineStyle='none')
xlabel('x')
ylabel('u')

subplot(3,1,3)
cla
plot(x, p, Color=[0.9290 0 0.1250], Marker='o', MarkerSize=3, LineStyle='none')
xlabel('x')
ylabel('p')
