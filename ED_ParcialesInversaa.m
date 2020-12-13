% RESOLUCION DE FRACCIONES PARCIALES E INVERSA DE LAPLACE
% ENTRADAS
    % Yf     -> Respuesta a la que se le quiere aplicar la inversa de laplace
% SALIDAS
    % Yf_inv -> Inversa de laplace de la respuesta
% Alejandra Soto y Susana Tristan    

function [Yf_inv] = ED_ParcialesInversaa(Yf)
syms t s Ys yo dyo Y c w; % Variables a utilizar
L = char(65:90);        % Variables a utilizar por las fracciones parciales
[Nf,Df] = numden((Yf)); % Numerador y denominador de la respuesta forzada

grado_Nf = polynomialDegree(Nf,s); % Grado del numerador
grado_Df = polynomialDegree(Df,s); % Grado del denominador

% Determinar si el grado del denominador es mayor al del numerador
if grado_Df > grado_Nf 
    DF = children(c*Df); % Grado del denominador mayor
    Bandera1 = 0;        % Aplicar resolucion por fracciones parciales
else
    [r,d] = polynomialReduce(Nf,Df,s); % Grado del numerador mayor o igual al del denominador  
    Yf = d + r/Df;                     % Realizar reducción del polinomio 
    Bandera1 = 1;                      % A la reduccion aplicar directo la inversa de laplace
end

if length(coeffs(Df,'all')) < 3 && length(children(numden(Yf))) ==1
   Bandera1=1; 
end    

if Bandera1 == 0
%%  OBTENER FRACCIONES PARCIALES
    cont_letra = 1;   % Contador de la variables 'letra' 
    DFF = factor(Df);
    Fracciones_f = 0;
    DFF_u = unique(DFF); % Denominadores unicos
    Bandera = 0;
    for i=1:length(DFF_u)
       lineal = isreal(roots(sym2poly(DFF(i)))); % Determina si es real o imaginarop (1 = real) 
       k=0;
       if lineal == 1 % lineal 
           for j=1:length(DFF)
               if DFF(j)== DFF_u(i)
                   k = k+1;
               end    
           end

           for j=1:k
              Fracciones_f = Fracciones_f + L(cont_letra)/(DFF_u(i)^j);
              cont_letra = cont_letra + 1;
           end 
    % Procedimiento para reales
    Bandera = 0;
       else
        % Imaginarios
         k=0;
           for j=1:length(DFF)
               if DFF(j)== DFF_u(i)
                   k = k+1;
               end    
           end

           for j=1:k
              Fracciones_f = Fracciones_f + (L(cont_letra)*s + L(cont_letra+1))/(DFF_u(i)^j);
              cont_letra = cont_letra + 1;
           end 
    % Procemientos para imaginarios
    Bandera =1;
       end    
    end 

    %% SEGUNDO DISPLAY 
    fprintf('RESOLUCION POR FRACCIONES PARCIALES \n')
    disp(Yf == Fracciones_f)

    %% RESOLVER FRACCIONES PARCIALES
    %% IMAGINARIOS
    if Bandera == 1
    Sol_1 = sym([]); % Division
    Sol_2 = sym([]); % Multiplicacion
    Sol_f = sym([]);
    disp_1 = 0; 
    disp_2 = 0;
    [Nf,Df] = numden((Yf)); % Numerador y denominador de la respuesta forzada
    Ff = children(expand(Fracciones_f));
    for i=1:length(Ff)
        [N_Frac,D_Frac] = numden(Ff(i));
        Sol_1(i) = 1/simplify(D_Frac/Df);
        disp_1 = disp_1 + N_Frac * Sol_1(i);
        Sol_2(i) = N_Frac*expand(Sol_1(i)); % Suma de Fracciones
        disp_2 = disp_2 + Sol_2(i);
    end
    Ff_collect = collect(disp_2);
    Ff_expand = expand(Sol_1);
    grados = [];
    for i = 1:length(Ff_expand)
        grados(i) = polynomialDegree(expand(Sol_2(i)),s);
    end
    m = max(grados);

    A = zeros(length(Sol_2),m+1);
    B = zeros(1,m+1);
    Ale=m+2;
    j=0;
    for i=1:length(Sol_2) % recorre renglones
        cn = s^(m+1) +  Nf;
        susy = polynomialDegree(Sol_2(i),s);
        if susy == 0
            cd = s^(m+1) +  Ff_expand(i);
        else 
            cd =  s^(m+1) + s^(susy)*Ff_expand(i);
        end
        coefsd = coeffs(cd,'All');
        coefsn = coeffs(cn,'All');
        A(i,:) = coefsd(2:end);
        B(1,i) = coefsn(Ale-j);
        j=j+1;
    end
    At=A';
    bt=B';
    Solucion = inv(At)*B';
    miguel = fliplr(Solucion');
    Solucion_final = 0 ;
    Ale=length(Solucion);
    for i=1:length(Solucion)
        [n,d]=numden(Ff(i));
          susy = polynomialDegree(n,s);
        Solucion_final =Solucion_final + (miguel(i)*s^susy)/d; 
    end  

    %% TERCER DISPLAY 
    disp(Nf == disp_1)
    disp(Nf == Ff_collect)
    m1_soluciones = children(Ff_collect);
    for i=1:length(m1_soluciones)
       m2_soluciones = children(c*m1_soluciones(i)); 
    end    
    fprintf('\n')
    for i=1:length(Solucion)
       fprintf(L(i))
       fprintf(' = ')
       fprintf(num2str(Solucion(i)))
       fprintf('\n')  
    end
    fprintf('\n')
    fprintf('Y(s) = ')
    fprintf(char(Solucion_final))
    fprintf('\n')
    fprintf('\n')
    else
    % REALES    
    Sol_1 = sym([]); % Division
    Sol_2 = sym([]); % Multiplicacion
    Sol_f = sym([]);
    disp_1 = 0; 
    disp_2 = 0;
    Ff = children(Fracciones_f);
    for i=1:length(Ff)
        [N_Frac,D_Frac] = numden(Ff(i));
        Sol_1(i) = 1/simplify(D_Frac/Df);
        disp_1 = disp_1 + N_Frac * Sol_1(i);
        Sol_2(i) = N_Frac*expand(Sol_1(i)); % Suma de Fracciones
        disp_2 = disp_2 + Sol_2(i);
    end
    Ff_collect = collect(disp_2);

    Ff_expand = expand(Sol_1);
    grados = [];
     for i = 1:length(Ff_expand)
        grados(i) = polynomialDegree(Ff_expand(i));
    end
    m = max(grados);
    A = zeros(length(Sol_2),m+1);
    B = zeros(1,m+1);
    for i=1:length(Sol_2) % recorre renglones
        cd = s^(m+1) +  Ff_expand(i);
        cn = s^(m+1) +  Nf;
        coefsd = coeffs(cd,'All');
        coefsn = coeffs(cn,'All');
        A(i,:) = coefsd(2:end);
        B(1,i) = coefsn(i+1);      
    end
    At=A';
    bt=B';
    Solucion = inv(At)*B';
    Solucion_final = 0 ;
    for i=1:length(Solucion)
        [n,d]=numden(Ff(i));
        Solucion_final =Solucion_final + (Solucion(i))/d; 
    end    

    %% TERCER DISPLAY 
    disp(Nf == disp_1)
    disp(Nf == Ff_collect)
    m1_soluciones = children(Ff_collect);
    for i=1:length(m1_soluciones)
       g = polynomialDegree(m1_soluciones(i),s); 
       if g > 0
       m2_soluciones = children(c*m1_soluciones(i)); 
       fprintf(char(m2_soluciones(end)))
       fprintf(' = ')
       fprintf(num2str(B(i)))
       fprintf('\n')
       else
       fprintf(char(m1_soluciones(end)))
       fprintf(' = ')
       fprintf(num2str(B(i)))
       fprintf('\n') 
       end    
    end    
    fprintf('\n')
    for i=1:length(Solucion)
       fprintf(L(i))
       fprintf(' = ')
       fprintf(num2str(Solucion(i)))
       fprintf('\n')  
    end
    fprintf('\n')
    fprintf('Y(s) = ')
    fprintf(char(Solucion_final))
    fprintf('\n')
    fprintf('\n')
    end
    %% INVERSA DE LAPLACE
    Yf_inv = ilaplace(Solucion_final);
    %% CUARTO DISPLAY 
    fprintf('INVERSA DE LAPLACE \n')
    fprintf('Y(t) = ')
    fprintf(char(Yf_inv))
    fprintf('\n')
else
    Yf_inv = ilaplace(Yf);
    %% CUARTO DISPLAY 
    fprintf('INVERSA DE LAPLACE \n')
    fprintf('Y(t) = ')
    fprintf(char(Yf_inv))
    fprintf('\n')
end  
end