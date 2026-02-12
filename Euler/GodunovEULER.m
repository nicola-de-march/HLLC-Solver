%% Godunov finite volume scheme for the shallow water equations (SWE) 
clear all 
close all
clc

global g fig1 fig2 gamma
global pdetype rm_solver

g = 9.81;                   % gravity constant on Earth 
gamma = 1.4;                % 
% Problems:
% 0: Shallow waters, 1: Euler
pdetype = 1;
% Type of Riemann solver
% 0: Exact, 1: HLLC, 2: Rusanov
rm_solver = 0;

%% Define the Riemann problem
hL   = 0.15;                % left water depth
hR   = 0.05;                % right water depth
uL   = 0;                   % left velocity
uR   = 0;                   % right velocity
psiL = 1;                   % left scalar 
psiR = 0;                   % right scalar 

QL = [hL; hL*uL; hL*psiL];  % left state vector
QR = [hR; hR*uR; hR*psiR];  % right state vector

time = 0;                   % initial and current time
tend = 0.5;                 % final time
xL = -1;                    % computational domain
xR = 1;  

IMAX = 100;                 % number of control volumes
dx = (xR-xL)/IMAX;          % mesh spacing 
x = linspace(xL+dx/2,xR-dx/2,IMAX);
CFL = 0.7;                  % Courant number 
NMAX = 10000;               % maximum number of time steps 

%% Define the initial condition
Q = zeros(3,IMAX);
for i=1:IMAX
    if (x(i)<=0) 
        Q(:,i) = QL;
    else
        Q(:,i) = QR;
    end
end
%% Plot initial condition 
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
title(sprintf('Time = %f',time))
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
%% Compute approximate solution using the finite volume scheme 
for n = 1:NMAX % main time loop
    % compute the maximum of the absolute value of all eigenvalues 
    amax = 0;
    for i=1:IMAX
        L = Lambda(Q(:,i));
        amax = max(amax,max(abs(L)));
    end
    % compute the time step using the CFL condition 
    dt = CFL*dx/amax;    
    % ensure that we reach tend
    if (time+dt>tend)
        dt = tend - time;
    end
    % stop criterion
    if (time>=tend)
        break
    end
    % loop over all spatial control volumes 
    Qnew = zeros(3,IMAX);
    for i=1:IMAX
        switch rm_solver
            case 0
                % Exact Riemann Solver
                if (i==1)
                    % Dirichlet BC on the left 
                    % Godunov flux
                    QGod = ExactRiemannSW(Q(:,i),Q(:,i+1),0);
                    fp = f(QGod);
                    QGod = ExactRiemannSW(2*QL-Q(:,i),Q(:,i),0);
                    fm = f(QGod);
                elseif (i==IMAX)
                     % Dirichlet BC on the right 
                    % Godunov flux
                    QGod = ExactRiemannSW(Q(:,i),2*QR-Q(:,i),0);
                    fp = f(QGod);
                    QGod = ExactRiemannSW(Q(:,i-1),Q(:,i),0);
                    fm = f(QGod);
                else
                    % Godunov flux
                    QGod = ExactRiemannSW(Q(:,i),Q(:,i+1),0);
                    fp = f(QGod);
                    QGod = ExactRiemannSW(Q(:,i-1),Q(:,i),0);
                    fm = f(QGod);
                end
            case 1
                % HLLC 
                if (i==1)
                    % Dirichlet BC on the left 
                    fp = HLLCSolver(Q(:,i),Q(:,i+1));
                    fm = HLLCSolver(2*QL-Q(:,i),Q(:,i));
                elseif (i==IMAX)
                     % Dirichlet BC on the right 
                    fp = HLLCSolver(Q(:,i),2*QR-Q(:,i));
                    fm = HLLCSolver(Q(:,i-1),Q(:,i));
                else
                    fp = HLLCSolver(Q(:,i),Q(:,i+1));
                    fm = HLLCSolver(Q(:,i-1),Q(:,i));
                end
            case 2
                % Rusanov
                if (i==1)
                    % Dirichlet BC on the left 
                    fp = Rusanov(Q(:,i),Q(:,i+1));
                    fm = Rusanov(2*QL-Q(:,i),Q(:,i));
                elseif (i==IMAX)
                     % Dirichlet BC on the right 
                    fp = Rusanov(Q(:,i),2*QR-Q(:,i));
                    fm = Rusanov(Q(:,i-1),Q(:,i));
                else
                    fp = Rusanov(Q(:,i),Q(:,i+1));
                    fm = Rusanov(Q(:,i-1),Q(:,i));
                end
        end
        % FV method
        Qnew(:,i) = Q(:,i) - dt/dx*(fp-fm);
    end  
    % Update time
    time = time + dt;
    % Overrite solution
    Q = Qnew;    
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
    title(sprintf('Time = %f',time))
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
    %% show the x-t diagram of the exact solution of the RP in a second figure 
    ierr = xtdiag(QL,QR,xL,xR,tend,time);
    drawnow
end

%% Plot exact solution, left part in blue, right part in yellow, as in the experiment 
fig1=figure(1);
% Compute solution with the exact Riemann solver on a finer mesh 
EMAX = 10*IMAX; 
xe = linspace(xL,xR,EMAX); 
Qe = zeros(3,EMAX);
for i= 1:EMAX
    xi = xe(i)/time;
    Qe(:,i) = ExactRiemannSW(QL,QR,xi);
end
psie   = Qe(3,:)./Qe(1,:); 
ue     = Qe(2,:)./Qe(1,:); 
lefte  = psie>0.5;
righte = psie<=0.5;
nl     = 0; 
for i=1:EMAX
    if(lefte(i)==1)
        nl = nl + 1; 
    end
end
subplot(3,1,1)
hold on
plot(xe(lefte), Qe(1,lefte),'b-')              
plot(xe(righte),Qe(1,righte),Color=[0.9290 0.6940 0.1250],LineStyle='-')
subplot(3,1,2)
hold on
plot(xe(lefte), ue(1,lefte),'b-')              
plot(xe(righte),ue(1,righte),Color=[0.9290 0.6940 0.1250],LineStyle='-')
subplot(3,1,3)
hold on
plot([xe(lefte),xe(nl)], [psie(1,lefte),0.5],'b-')              
plot([xe(nl),xe(righte)],[0.5,psie(1,righte)],Color=[0.9290 0.6940 0.1250],LineStyle='-')
drawnow