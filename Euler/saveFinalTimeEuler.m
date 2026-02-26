%% Save final time
function saveFinalTimeEuler(Q,x,t, exact)

global rm_solver

% Converte Q in variabili primitive
V = PDECons2Prim(Q);  % V = [rho; u; p]

rho = V(1,:);
u   = V(2,:);
p   = V(3,:);

% 0: Exact, 1: HLLC, 2: Rusanov, 3: TV Splitting
if (exact)
    rm_solver_string = 'Exact';
else
    switch rm_solver    
        case 0
            rm_solver_string = 'Exact';
        case 1
            rm_solver_string = 'HLLC';
        case 2
            rm_solver_string = 'Rusanov';
        case 3
            rm_solver_string = 'TV';
        otherwise
            error('Unknown Riemann solver.')
    end
end
% Create data structure
data.x   = x;
data.t   = t;
data.rho = rho;
data.u   = u;
data.p   = p;
data.solver = rm_solver_string;

% Create folder if it does not exist
folder_name = fullfile('Output', rm_solver_string);
if ~exist(folder_name, 'dir')
    mkdir(folder_name);
end

% File name (time with fixed formatting to avoid floating errors)
file_name = sprintf('Euler_%s.mat', rm_solver_string);
full_path = fullfile(folder_name, file_name);

% Save
save(full_path, 'data');