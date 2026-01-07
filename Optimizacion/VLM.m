function [C_am,Emax,CLopt,flecha_14_rad,flecha_12_rad,t_r] = VLM(S,c,lambda,A,flecha,CD0,grafica)
    %% PERFIL NACA 3414
    fmax = 0.03;                       
    xfmax = 0.4;                              
    
    %Pendiente del perfil para 3c/4
    dZc2 = @(x) (2*fmax/((1-xfmax)^2)*(xfmax-x));
    dzdx = dZc2(3*c/4);
    
    %% DATOS (ANEXO 1)
    epsilon_grados = 0;
    epsilon = deg2rad(epsilon_grados);
    Ny = 20;
    Nx = 1;
    y1 = 0; 
    alpha_grados = linspace(-5,25,100)';
    alpha = deg2rad(alpha_grados);
    s = size(alpha,1);
    U_inf = 1;
    
    %% CUERDAS Y PARAMETROS (ANEXO 2)
    b = sqrt(A*S);
    C_raiz = S / ((b/2)*(1+lambda)+y1*(1-lambda));
    C_punta = lambda*C_raiz;
    flecha_ba_rad = (atan(tan(flecha)+(C_raiz-C_punta)/(4*((b/2)-y1))));
    flecha_bs_rad = (atan((C_raiz-C_punta)/((b/2)-y1)-tan(flecha_ba_rad)));
    C_am = (2/3)*C_raiz*(1+lambda+lambda^2)/(1+lambda);

    %% COORDENADAS 8 PUNTOS (ANEXO 3)
    ala = zeros(2,8); %filas = coordenada x,y 
                      %columnas = puntos
    ala(1,1) = b/2*tan(flecha_ba_rad);
    ala(2,1) = b/2;
    ala(1,2) = 0;
    ala(2,2) = 0;
    ala(1,3) = b/2*tan(flecha_ba_rad) + C_punta;
    ala(2,3) = b/2;
    ala(1,4) = C_raiz;
    ala(2,4) = 0;
    ala(1,5) = b/2*tan(flecha_ba_rad) + 1/4*C_punta;
    ala(2,5) = b/2;
    ala(1,6) = 1/4*C_raiz;
    ala(2,6) = 0;
    ala(1,7) = b/2*tan(flecha_ba_rad)+3/4*C_punta;
    ala(2,7) = b/2;
    ala(1,8) = 3/4*C_raiz;
    ala(2,8) = 0;
    
    %% ANGULOS (ANEXO 4)
    flecha_14_rad = atan((b/2*tan(flecha_ba_rad) + 1/4*C_punta - 1/4*C_raiz) / (b/2));
    flecha_34_rad = atan((b/2*tan(flecha_ba_rad)+3/4*C_punta - 3/4*C_raiz) / (b/2));
    flecha_12_rad = atan((b/2*tan(flecha_ba_rad) + 1/2*C_punta - 1/2*C_raiz) / (b/2));
    t_r = 0.14*C_raiz;
    %% COORDENADAS ESQUINAS PANELES (ANEXO 5)
    %Coordenada x
    for j = 1:Ny+1 %1:21
        for i = 1:Nx+1 
            if i==1
                coord_paneles_dcha_x(i,j) = b/(2*Ny)*(j-1)*tan(flecha_ba_rad);
            end
            if i==2
                coord_paneles_dcha_x(i,j) = C_raiz - b/(2*Ny)*(j-1)*tan(flecha_bs_rad);
            end
        end
    end
    
    %Coordenada y
    for j = 1:Ny+1
        for i = 1:Nx+1
            coord_paneles_dcha_y(i,j) = b/(2*Ny)*(j-1);
        end
    end
    
    %Esquinas de los paneles DERECHA (coloco los vectores anteriores)
    paneles_dcha = zeros(Ny*2,4); %hay Ny paneles, con 4 esquinas cada uno, la  matriz esta formada por 2 filas para cada panel de forma que en la primera fila de cada panel esta la coordenada x de cada esquina y en la segunda la coordenada y, en total habra 4 columnas para cada esquina, nombradas desde la izq inferior, en sentido antihorario
    paneles_izq = zeros(Ny*2,4); 
    k=0;
    
    for i=1:Ny*2
        %coloco las coordenadas x en filas impares------------------
        if rem(i,2)~= 0 && i~=1  %contador para que ubique el panel mas abajo, p
            k=k+1;
        end
    
        for j=1:4
            if rem(i,2)~= 0  %si i es impar es que es fila de coordenada x, uso vector de coordenadas de x
                p=i-k;
                paneles_dcha(i,1) = coord_paneles_dcha_x(1,p);
                paneles_dcha(i,2) = coord_paneles_dcha_x(1,p+1);
                paneles_dcha(i,3) = coord_paneles_dcha_x(2,p);
                paneles_dcha(i,4) = coord_paneles_dcha_x(2,p+1);
    
                paneles_izq(i,1) = coord_paneles_dcha_x(1,p);
                paneles_izq(i,2) = coord_paneles_dcha_x(1,p+1);
                paneles_izq(i,3) = coord_paneles_dcha_x(2,p);
                paneles_izq(i,4) = coord_paneles_dcha_x(2,p+1);
                    
            end
        %coloco las coordenadas y en filas pares------------------
            if rem(i,2)==0 %%si i es par es que es fila de coordenada y, uso vector de coordenadas de y
                p=i-i/2;
                paneles_dcha(i,1) = coord_paneles_dcha_y(1,p);
                paneles_dcha(i,2) = coord_paneles_dcha_y(1,p+1);
                paneles_dcha(i,3) = coord_paneles_dcha_y(2,p);
                paneles_dcha(i,4) = coord_paneles_dcha_y(2,p+1); 
    
                paneles_izq(i,1) = -coord_paneles_dcha_y(1,p);
                paneles_izq(i,2) = -coord_paneles_dcha_y(1,p+1);
                paneles_izq(i,3) = -coord_paneles_dcha_y(2,p);
                paneles_izq(i,4) = -coord_paneles_dcha_y(2,p+1); 
            end
        end
    end
    
    %% COORDENADAS DE PUNTOS DE CONTROL (ANEXO 6)
    pc_dcha = zeros(2,Ny);
    pc_izq = zeros(2,Ny);
    
    for i = 1:2
        p = 0;
        for j = 1:Ny  %numero de paneles
            %coordenadas y
            p=p+1; %panel
            k=p-1;
            pc_dcha(2,j) = (b/(2*Ny))*1/2 *(p+k);
            pc_izq(2,j) = -(b/(2*Ny))*1/2 *(p+k);
            %coordenadas x
            pc_dcha(1,j) = C_raiz*3/4 + pc_dcha(2,j) * tan(flecha_34_rad);
            pc_izq(1,j) = C_raiz*3/4 - pc_izq(2,j) * tan(flecha_34_rad); %el - para que se haga postitivo como en dcha   
        end
    end
    
    %% COORDENADAS DE ESQUINAS TORBELLINOS (ANEXO 7)
    esquina_torb_dcha = zeros(Ny*2,2);
    esquina_torb_izq = zeros(Ny*2,2);
    
    %coordenadas en vectores separados
    coord_torb_dcha = zeros(2,Ny+1);
    coord_torb_izq = zeros(2,Ny+1);
    
    for j=1:Ny+1
        coord_torb_dcha(2,j) = b/(2*Ny)*(j-1);
        coord_torb_izq(2,j) = - b/(2*Ny)*(j-1);
        coord_torb_dcha(1,j) = C_raiz/4 + coord_torb_dcha(2,j)*tan(flecha_14_rad);
        coord_torb_izq(1,j) = C_raiz/4 - coord_torb_izq(2,j)*tan(flecha_14_rad);  %pongo - para que quede positivo
    end
    
    k=0;
    for i=1:Ny*2
        %coloco las coordenadas x en filas impares------------------
        if rem(i,2)~= 0 && i~=1  %contador para que ubique el panel mas abajo, p
            k=k+1;
        end
    
        for j=1:2
            if rem(i,2)~= 0  %si i es impar es que es fila de coordenada x, uso vector de coordenadas de x
                p=i-k;
                esquina_torb_dcha(i,1) = coord_torb_dcha(1,p);
                esquina_torb_dcha(i,2) = coord_torb_dcha(1,p+1);
    
                esquina_torb_izq(i,1) = coord_torb_izq(1,p);
                esquina_torb_izq(i,2) = coord_torb_izq(1,p+1);
            end
        %coloco las coordenadas y en filas pares------------------
            if rem(i,2)==0 %%si i es par es que es fila de coordenada y, uso vector de coordenadas de y
                p=i-i/2;
                esquina_torb_dcha(i,1) = coord_torb_dcha(2,p);
                esquina_torb_dcha(i,2) = coord_torb_dcha(2,p+1);
    
                esquina_torb_izq(i,1) = coord_torb_izq(2,p);
                esquina_torb_izq(i,2) = coord_torb_izq(2,p+1);
            end
        end
    end
    
    %% ENVERGADURA DE PANELES Y CUERDA MEDIA GEOMETRICA (ANEXO 8)
    envergadura_paneles = b/(2*Ny);  %en todos los paneles se tiene la misma envergadura ya que viene de una particion equidistante de la envergadura total del ala
    cgm_paneles_semiala = zeros(Ny,1);  %ambas semialas tienen las mismas cgm, solo se calculan los valores de una de ellas
    
    cuerda = zeros(Ny,1);
    m = 1;
    for i = 1:Ny
        cuerda(i,1) = abs((paneles_dcha(m,3)-paneles_dcha(m,1)));
        cgm_paneles_semiala(i,1) = (cuerda(i,1) + abs((paneles_dcha(m,4)-paneles_dcha(m,2))))/2;
        m = m + 2;
    end
    
    sj = envergadura_paneles*cgm_paneles_semiala;
    
    %% COORDENADAS DE 1/4 PANELES
    c14_dcha = zeros(2,Ny);
    c14_izq = zeros(2,Ny);
    
    for i = 1:2
    p=0;
        for j = 1:Ny  %numero de paneles
            %coordenadas y
            p=p+1; %panel
            k=p-1;
            c14_dcha(2,j) = (b/(2*Ny))*1/2 *(p+k);
            c14_izq(2,j) = -(b/(2*Ny))*1/2 *(p+k);
            %coordenadas x
            c14_dcha(1,j) = C_raiz*1/4 + c14_dcha(2,j) * tan(flecha_14_rad);
            c14_izq(1,j) = C_raiz*1/4 - c14_izq(2,j) * tan(flecha_14_rad); %el - para que se haga postitivo como en dcha
            
        end
    end
    
    %% CÁLCULO DE LOS TORBELLINOS
    % i: panel inducido, j: panel inductor
    a_dcha = zeros(Ny); a_izq = zeros(Ny);
    b_dcha = zeros(Ny); b_izq = zeros(Ny);
    c_dcha = zeros(Ny); c_izq = zeros(Ny);
    d_dcha = zeros(Ny); d_izq = zeros(Ny);
    e_dcha = zeros(Ny); e_izq = zeros(Ny);
    f_dcha = zeros(Ny); f_izq = zeros(Ny);
    g_dcha = zeros(Ny); g_izq = zeros(Ny);
    h_dcha = zeros(Ny); h_izq = zeros(Ny);
    k_dcha = zeros(Ny); k_izq = zeros(Ny);
    l_dcha = zeros(Ny); l_izq = zeros(Ny);
    V_dcha = zeros(Ny); V_izq = zeros(Ny);
    
    for i = 1:Ny
       m = 1;
       for j = 1:Ny
           if m <= 2*Ny
               %Ala derecha
               a_dcha(i,j) = pc_dcha(1,i) - esquina_torb_dcha(m,1);
               b_dcha(i,j) = pc_dcha(2,i) - esquina_torb_dcha(m+1,1);
               c_dcha(i,j) = pc_dcha(1,i) - esquina_torb_dcha(m,2);
               d_dcha(i,j) = pc_dcha(2,i) - esquina_torb_dcha(m+1,2);
               e_dcha(i,j) = sqrt(a_dcha(i,j)^2 + b_dcha(i,j)^2);
               f_dcha(i,j) = sqrt(c_dcha(i,j)^2 + d_dcha(i,j)^2);
               g_dcha(i,j) = esquina_torb_dcha(m,2) - esquina_torb_dcha(m,1);
               h_dcha(i,j) = esquina_torb_dcha(m+1,2) - esquina_torb_dcha(m+1,1);
               k_dcha(i,j) = (g_dcha(i,j)*a_dcha(i,j) + h_dcha(i,j)*b_dcha(i,j))/e_dcha(i,j) - (g_dcha(i,j)*c_dcha(i,j) + h_dcha(i,j)*d_dcha(i,j))/f_dcha(i,j);
               l_dcha(i,j) = -(1 + a_dcha(i,j)/e_dcha(i,j))/(b_dcha(i,j)) + (1/d_dcha(i,j))*(1 + c_dcha(i,j)/f_dcha(i,j));
               V_dcha(i,j) = (k_dcha(i,j)/(a_dcha(i,j)*d_dcha(i,j) - c_dcha(i,j)*b_dcha(i,j)) + l_dcha(i,j))/(4*pi);
               %Ala izquierda
               a_izq(i,j) = pc_dcha(1,i) - esquina_torb_izq(m,2);
               b_izq(i,j) = pc_dcha(2,i) - esquina_torb_izq(m+1,2);
               c_izq(i,j) = pc_dcha(1,i) - esquina_torb_izq(m,1);
               d_izq(i,j) = pc_dcha(2,i) - esquina_torb_izq(m+1,1);
               e_izq(i,j) = sqrt(a_izq(i,j)^2 + b_izq(i,j)^2);
               f_izq(i,j) = sqrt(c_izq(i,j)^2 + d_izq(i,j)^2);
               g_izq(i,j) = esquina_torb_izq(m,1) - esquina_torb_izq(m,2);
               h_izq(i,j) = esquina_torb_izq(m+1,1) - esquina_torb_izq(m+1,2);
               k_izq(i,j) = (g_izq(i,j)*a_izq(i,j) + h_izq(i,j)*b_izq(i,j))/e_izq(i,j) - (g_izq(i,j)*c_izq(i,j) + h_izq(i,j)*d_izq(i,j))/f_izq(i,j);
               l_izq(i,j) = -(1 + a_izq(i,j)/e_izq(i,j))/(b_izq(i,j)) + (1/d_izq(i,j))*(1 + c_izq(i,j)/f_izq(i,j));
               V_izq(i,j) = (k_izq(i,j)/(a_izq(i,j)*d_izq(i,j) - c_izq(i,j)*b_izq(i,j)) + l_izq(i,j))/(4*pi);
               m = m + 2;
          end
       end
    end
    
    a = zeros(Ny);
    bj = zeros(Ny,s);
    gamma = zeros(Ny,s);
    
    for n = 1:s
        for i = 1:Ny
            for j = 1:Ny
                a(i,j) = V_dcha(i,j) + V_izq(i,j);
            end
            bj(i,n) = alpha(n) + epsilon*(pc_dcha(2,i))/(b/2) - dzdx;
        end
        gamma(:,n) = a\(-bj(:,n));
    end
    
    %% COEFICIENTE DE SUSTENTACIÓN
    cl = zeros(Ny,s);
    CL = zeros(1,s);
    
    for n = 1:s
        sum = 0;
        for j = 1:Ny
            cl(j,n) = 2*gamma(j,n)/(cgm_paneles_semiala(j)*U_inf);
            sum = sum + cl(j,n)*envergadura_paneles*cgm_paneles_semiala(j);
        end
        CL(n) = (2/S)*sum;
    end
    
    cl_basica = zeros(Ny,1);
    cl_adicional = zeros(Ny,1);
    
    for j = 1:Ny
        cl_adicional(j,1) = (cl(j,1) - cl(j,2))/(CL(1) - CL(2));
        cl_basica(j,1) = cl(j,1) - cl_adicional(j,1)*CL(1);
    end
    
    %% COEFICIENTE DE RESISTENCIA INDUCIDA
    wlibre_dcha = zeros(Ny,Ny);
    wlibre_izq = zeros(Ny,Ny);
    gamma_T = zeros(Ny,s);
    alpha_inf = zeros(Ny,s);
    cd = zeros(Ny,s);
    CDi = zeros(1,s);
    
    for n = 1:s
        m = 1;
        for i = 1:Ny
            for j = 1:Ny
                if j < Ny
                    gamma_T(j,n) = gamma(j,n) - gamma(j+1,n);
                    wlibre_dcha(i,j) = (1/(2*pi))/(pc_dcha(2,j) - esquina_torb_dcha(m+1,2));
                    wlibre_izq(i,j) = -(1/(2*pi))/(pc_dcha(2,j) - esquina_torb_izq(m+1,2));
                else
                    gamma_T(j,n) = gamma(j,n);
                    wlibre_dcha(i,j) = (1/(2*pi))/(pc_dcha(2,j) - esquina_torb_dcha(m+1,2));
                    wlibre_izq(i,j) = -(1/(2*pi))/(pc_dcha(2,j) - esquina_torb_izq(m+1,2));
                end
            end
            m = m + 2;
            wlibre = wlibre_dcha + wlibre_izq;
            alpha_inf = (gamma_T'*wlibre/U_inf)';
        end
    end
    
    for n = 1:s
        sum = 0;
        for j = 1:Ny
            cd(j,n) = -cl(j,n)*alpha_inf(j,n)/2;
            sum = sum + cd(j,n)*envergadura_paneles*cgm_paneles_semiala(j);
        end
        CDi(n) = (2/S)*sum;
    end

    % Cálculo de la eficiencia máxima
    CD = CD0 +CDi;
    E = CL./CD;
    [Emax, idx] = max(E);
    alphaeff = alpha_grados(idx);
    CLopt = CL(idx);
    CDopt = CD(idx);
        
    %% Figuras
    y_ext_dcha = [paneles_dcha(2,1:2:end), fliplr(paneles_dcha(4,1:2:end))];
    x_ext_dcha = [paneles_dcha(1,1:2:end), fliplr(paneles_dcha(3,1:2:end))];

    % --- Ala izquierda ---
    y_ext_izq = [paneles_izq(2,1:2:end), fliplr(paneles_izq(4,1:2:end))];
    x_ext_izq = [paneles_izq(1,1:2:end), fliplr(paneles_izq(3,1:2:end))];

    y_ext_dcha = [y_ext_dcha, y_ext_dcha(1)];
    x_ext_dcha = [x_ext_dcha, x_ext_dcha(1)];
    
    % Ala izquierda
    y_ext_izq = [y_ext_izq, y_ext_izq(1)];
    x_ext_izq = [x_ext_izq, x_ext_izq(1)];
    
    % %Gráfica panelización
    if grafica
        figure;
        plot([y_ext_dcha, y_ext_dcha(1), y_ext_izq, y_ext_izq(1)], ...
        [x_ext_dcha, x_ext_dcha(1), x_ext_izq, x_ext_izq(1)]);
        xlabel('y [m]');                                                      
        ylabel('x [m]');
        title('Contorno exterior del ala');
        hold off
        
        %Gráfica CL vs. alpha
        figure;
        plot(alpha_grados,CL,'color',[0.91,0.46,0.26])
        xlabel('\alpha [º]');                                                      
        ylabel('C_L [-]');
        title('Coeficiente de sustentación en función de \alpha');
    
        %Gráfica CD vs. alpha
        CD = CD0+CDi;
        figure;
        plot(alpha_grados,CD,'color',[0.48,0.71,0.32])
        xlabel('\alpha [º]');                                                      
        ylabel('C_D [-]');
        title('Coeficiente de resistencia en función de \alpha');

        % Gráfica CL/CD vs. alpha
        figure;
        plot(alpha_grados,E,'color',[0.48,0.71,0.32])
        xlabel('\alpha [º]');                                                      
        ylabel('E [-]');
        title('Eficiencia aerodinámica en función de \alpha');
        
        fprintf('Aeródinamica\n');
        fprintf('Emax = %.4f para alpha = %.4f\n', Emax, alphaeff);
        fprintf('CL = %.4f\n', CLopt);
        fprintf('CD = %.4f\n', CDopt);
    end
end
