% =========================================================================
% SOLVER FOR 1D SHALLOW WATER AND EULER'S EQUATION
% =========================================================================
clear all
close all
clc
%% Add path
addpath('Euler');
addpath('ShallowWaters');
addpath('FluxSolvers');
addpath('Utils');

%% Settup problem and solver

global pdetype rm_solver 

% PROBLEM TO SOLVE
% 0: Shallow waters, 1: Euler
pdetype = 1;
% TYPE OF RIEMANN SOLVERS
% 0: Exact, 1: HLLC, 2: Rusanov, 3: TV Splitting
rm_solver = 3;

% Define global parameters
global g gamma
g = 9.81; 
gamma = 1.4;

%% Init grids and time steps
time = 0;                   % initial and current time
tend = 0.5;                 % final time
xL = -1;                    % computational domain
xR = 1;  

IMAX = 100;                 % number of control volumes
dx = (xR-xL)/IMAX;          % mesh spacing 
x = linspace(xL+dx/2,xR-dx/2,IMAX);
CFL = 0.9;                  % Courant number 
NMAX = 10000;               % maximum number of time steps 

%% Initial condition: init solution Q
switch pdetype  
    case 0
        InitSWE;
        plotSWE(Q, x, 0)
    case 1
        InitEuler;
        plotEuler(Q, x, 0)
end

%% Solve with Finite Voleume scheme

for n = 1:NMAX
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

    % Solve
    switch pdetype
        case 0
            GodunovSWE;
        case 1
            GodunovEULER;
    end
    % Update time
    time = time + dt;
    % Overrite solution
    Q = Qnew;  
    switch pdetype  
        case 0
            plotSWE(Q, x, time);
            drawnow
            hold on
        case 1
            plotEuler(Q, x, time);
            drawnow
    end
end

%% Plot exact solution
switch pdetype
    case 0
        plotExactSWE;
    case 1
        plotExactEuler;
end






