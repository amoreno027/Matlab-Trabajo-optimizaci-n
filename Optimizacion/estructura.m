function [Wi,Wf,V] = estructura(lambda,b,t_c,S,A,g,rho,CLopt,rho_fuel,rho_material,K_ala,eta_ult,W0,flecha_14_rad,flecha_12_rad,c_am,t_r)
    %fprintf('\nEstructura\n');
    %fprintf('Peso vacío operativo sin alas: W0 = %.4f\n', W0); 

    % Peso del combustible
    tau = 1;
    Vfuel = 0.54*S^2/b*t_c*(1+lambda*sqrt(tau)+lambda^2*tau)/(1+lambda)^2; % Volumen disponible de combustible
    Wfuel = g*rho_fuel*Vfuel;
    %fprintf('Peso del combustible: Wf = %.4f\n', Wfuel);

    % Peso estructural del ala
    Ww = S*c_am*t_c*rho_material*K_ala*(A*eta_ult/cos(flecha_14_rad))^0.6*lambda^0.04*g;
    %fprintf('Peso de las alas: Ww = %.4f\n', Ww);
    
    % Peso total
    Wi = W0 + Ww + Wfuel; 
    %fprintf('Peso total inicial: Wi = %.4f\n', Wi);

    % Velocidad crucero
    Wf = Wi - Wfuel;
    Wmedio = 1/2*(Wf + Wi);
    V = sqrt(Wmedio/(1/2*rho*S*CLopt));
    %fprintf('Velocidad crucero: V = %.4f\n', V);