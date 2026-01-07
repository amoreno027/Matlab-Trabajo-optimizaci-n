close all; clc;
    fprintf('Optimización del alcance del UAV\n');

    % Punto inicial
    x0 = [1, sqrt(13), 0.15];   % [S,b,t/c]
    lb = [1; 2; 0.1];       % Límite inferior
    ub = [3; 6; 0.25];      % Límite superior

    % Modo análisis
    R = -fobj(x0,true); % True = plots, False = no plots
    fprintf('Alcance: %.1f km\n',R)

    % Modo optimización
    % Método del gradiente
    algos = {'sqp', 'interior-point', 'active-set'}; % Para compararlos y ver cuál utilizar
    results = struct();  % estructura para guardar resultados

    for k = 1:length(algos)
        options = optimoptions('fmincon','Display','iter','Algorithm',algos{k});
        [x_opt, fval, exitflag, output] = fmincon(@(x) fobj(x, false), x0, [], [], [], [], lb, ub, [], options);
    
        % Guardar resultados
        results(k).algorithm = algos{k};
        results(k).x_opt = x_opt;
        results(k).R_max = -fval;
        results(k).iterations = output.iterations;
    
        fprintf('Algoritmo: %s\n', algos{k});
        fprintf(' Diseño óptimo: [S b t/c] = %s\n', mat2str(x_opt,4));
        fprintf(' Alcance máximo: %.1f km\n', -fval);
        fprintf(' Iteraciones: %d\n\n', output.iterations);
        
        fobj(x_opt,true);
    end

    % Método heurístico

    % Optimización multiobjetivo



