function f = fobj(x,grafica)
    %% VARIABLES DE DISEÑO
    S = x(1);
    b = x(2);
    t_c = x(3);
    
    %% PARÁMETROS FIJOS
    g = 9.81; % [m/s^2]
    rho = 1.225; % [kg/m^3]
    
    CD0 = 0.03; % Coeficiente de resistencia parásita
    lambda = 0.5; % Estrechamiento
    
    m0 = 20; % Masa en vacío (sin el ala) [kg]
    
    cT = 1.5/3600; % Consumo específico del motor [s^-1]
    rho_fuel = 800; % Densidad del combustible [kg/m^3]
    
    rho_material = 2780; % Densidad del material de las alas (Aluminio aeronáutico (2024 / 7075)) [kg/m^3]
    K_ala = 0.00125; % Factor de densidad del ala
    eta_ult = 1.5*1.75; % Factor de carga último
    
    flecha = 0; % [º]
    flecha_rad = deg2rad(flecha);
    
    %% RELACIONES DEL MODELO
    A = b^2/S; % Alargamiento
    c = S/b; % Cuerda del ala
    W0 = m0*g; % Peso vacío operativo sin ala
    
    %% CÁLCULO DE AERODINÁMICA
    [c_am,Emax,CLopt,flecha_14_rad,flecha_12_rad,t_r] = VLM(S,c,lambda,A,flecha_rad,CD0,grafica);
    [Wi,Wf,V] = estructura(lambda,b,t_c,S,A,g,rho,CLopt,rho_fuel,rho_material,K_ala,eta_ult,W0,flecha_14_rad,flecha_12_rad,c_am,t_r);
    R = alcance(Emax,V,Wi,Wf,cT);
    f = -R;
end