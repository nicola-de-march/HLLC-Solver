clear
close all
clc

results_folder = 'Output';

%% Detect methods
dir_info = dir(results_folder);
dir_info = dir_info([dir_info.isdir]);
dir_info = dir_info(~ismember({dir_info.name}, {'.','..'}));

methods = string({dir_info.name});

fprintf('Detected methods:\n')
disp(methods)


%% Filter methods
% 
tmp = ["Rusanov","Exact","TV"];
methods = methods(ismember(methods, tmp));

%% Plot
figure('Position',[300 200 800 800])
tiledlayout(3,1)

for m = 1:length(methods)

    method = methods(m);

    folder_name = fullfile(results_folder, method);
    file_name   = sprintf('Euler_%s.mat', method);
    full_path   = fullfile(folder_name, file_name);

    if ~exist(full_path, 'file')
        warning('File not found: %s', full_path);
        continue
    end

    load(full_path);   % loads struct "data"


    % Choose style depending on method
    if method == "Exact"
        lineStyle = '-';
        marker    = 'none';
        lineWidth = 2.5;
    else
        lineStyle = 'none';
        lineWidth = 1;
        
        
        % Different markers for different methods
        switch method
            case "HLLC"
                marker = 'o';
            case "Rusanov"
                marker = 'x';
            case "TV"
                marker = 's';
            otherwise
                marker = 'd';
        end
    end

    markersize = 6;
    
    % Density
    nexttile(1)
    hold on
    plot(data.x, data.rho, ...
        'LineStyle', lineStyle, ...
        'Marker', marker, ...
        'MarkerSize', markersize, ...
        'LineWidth', lineWidth)
    
    % Velocity
    nexttile(2)
    hold on
    plot(data.x, data.u, ...
        'LineStyle', lineStyle, ...
        'Marker', marker, ...
        'MarkerSize', markersize, ...
        'LineWidth', lineWidth)
    
    % Pressure
    nexttile(3)
    hold on
    plot(data.x, data.p, ...
        'LineStyle', lineStyle, ...
        'Marker', marker, ...
        'MarkerSize', markersize, ...
        'LineWidth', lineWidth)

end

nexttile(1)
title(sprintf('Final time = %.6f\n', data.t))
ylabel('\rho')
legend(methods, 'Location','best')
grid on

nexttile(2)
ylabel('u')
grid on

nexttile(3)
ylabel('p')
xlabel('x')
grid on