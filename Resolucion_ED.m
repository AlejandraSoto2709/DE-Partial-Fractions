% RESOLUCION DE ECUACIONES DIFERENCIALES CON CONDICIONES INICIALES
% LAPLACE Y FRACCIONES PARCIALES 
% NOTAS: 
    % ED -> d^2/dt^2 y(t) + d/dt y(t) + y(t) = d/dt x(t) + x(t) (Se insertan como polinomios)
    % L{y}   -> Y(s)                   (Representado Ys)
    % L{y'}  -> sY(s) - y(0)           (Representado: sYs - yo) 
    % L{y''} -> s^2Y(s) - sy(0) -y'(0) (Representado: s^2Ys - syo - dyo)
 % ENTRADAS:
    % p1  -> Parte 1 de la ecuacion diferencial (Antes del igual) 
    % p2  -> Parte 2 de la ecuacion diferencial (Despues del igual) 
    % y0  -> (Valor en y = 0)
    % x_t -> x(t) 
    % dy0 -> Valor de la derivada en y = 0
 % SALIDAS:
    %  EC_3    -> Despeje de Y(s)
    %  Yf      -> Respuesta forzada Yf(s)
    %  Yf_inv  -> Inversa de laplace de la respuesta forzada Yf(t)
    %  Yn      -> Respuesta natural Yn(s) 
    %  Yn_inv  >  Inversa de laplace de la respuesta natural Yn(t)
% Alejandra Soto y Susana Tristan 

function [EC_3,Yf,Yn,Yf_inv,Yn_inv] = Resolucion_ED(p1,p2,y0,x_t,dy0)
syms t s Ys yo dyo Y c w % Varaibles simbolicas a utilizar
L = char(65:90);         % Variables que se utilizar en fracciones parciales 

%% LAPLACE y, y', y''
L_y = Ys;
L_d1y = s*Ys - yo;
L_d2y = s^2*Ys - s*yo - dyo;

%% LAPLACE
P1 = sym(p1); % Primera parte de la ED. en simbolico
P2= sym(p2); % Segunda parte de la ED. en simbolico

% Laplace de la primera parte.
L_P1 = simplify(P1(1).*L_d2y+ P1(2).*L_d1y + P1(3).*L_y);
% Laplace de la segunda parte.
L_P2 = simplify(P2(1)*s*laplace(x_t)+P2(2)*laplace(x_t));
% Sustitucion de condiciones iniciales 
L_P1_cond = subs(L_P1,{yo,dyo},{y0,dy0});
% Factorizacion de Ys.
EC = L_P1_cond == L_P2; % Ec. con Laplace y condiciones iniciales
EC_1 = collect(EC,Ys);  % Agrupacion de terminos con Ys
Ch = children(EC_1);
Ch_2 = children(Ch(1));

n = 1; p = 0;
while n < length(Ch_2)
   p = p - Ch_2(n + 1);
   n = n + 1; 
end

EC_2 = Ch_2(1) == Ch(2) + p;
Ch_3 = children(Ch_2(1));
Ch_4 = children(EC_2);
EC_3 = Ch_3(2) ==  Ch_4(2)/Ch_3(1); % Despeje Ys

%%   RESPUESTA NATURAL Y FORZADA
Ch_5 = children(Ch_2(1));
A = Ch_5(1); 
B = Ch_2(2);
C = Ch(2);
% Respuesta natural 
Yn = (-1*B)/A;
% Respuesta forzada 
Yf = C/A;

%% PRIMER DISPLAY  
disp('LAPLACE')
disp(L_P1 == L_P2)
disp(EC)
disp(EC_1)
disp(EC_2)
disp(EC_3)
fprintf( '\n')
fprintf('Respuesta \n')
disp(Ys == Yf + Yn)
fprintf('Respuesta forzada \n')
fprintf('Yf = ' )
disp(Yf)
fprintf('Respuesta natural \n')
fprintf('Yn = ')
disp(Yn)

%% SOLUCION POR FRACCIONES PARCIALES E INVERSA DE LAPLACE
% Respuesta Forzada
disp('RESPUESTA FORZADA')
Yf_inv = ED_ParcialesInversaa(Yf);
fprintf('\n')
fprintf('\n')

% Respuesta Natural 
disp('RESPUESTA NATURAL')
Yn_inv = ED_ParcialesInversaa(Yn);
fprintf('\n')
fprintf('\n')

% Respuesta Completa
Yc = Yf_inv + Yn_inv;
disp('RESPUESTA')
fprintf('Y(t) = ')
disp(Yc)
fprintf('\n')
