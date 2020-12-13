%%  RESOLUCION DE ECUACIONES DIFERENCIALES CON CONDICIONES INICIALES
% LAPLACE Y FRACCIONES PARCIALES 
% NOTAS : 
    % ED -> d^2/dt^2 y(t) + d/dt y(t) + y(t) = d/dt x(t) + x(t) (Se insertan como polinomios)
    % Encuentra la inversa de laplace de la version extendida de las fracciones parciales 
% Alejandra Soto y Susana Tristan

%% DATOS DE ED. 
clear all, close all, clc
syms t s Ys yo dyo Y c w % Varaibles simbolicas a utilizar
L = char(65:90);         % Variables que se utilizar en fracciones parciales 

% Señales y Sistemas - Haykin
% Ejercicios 6.38 
% a)
    p1 = [0,1,10]; % Primera parte de la ED.
    p2 = [0,10];        % Segunda parte de la ED. [p1 = p2]
    y0 = 1;             % Condicion inical 
    x_t = heaviside(t); % Valor de x(t)
    dy0 = 0;
% b)
%    p1 = [1,5,6];
%    p2 = [-3,-4];
%    y0 = -1;
%    x_t = exp(-t) * heaviside(t);
%    dy0=5;
%  c)
%     p1 = [1,1,0];
%     p2 = [0,8];
%     y0 = 0;
%     x_t = exp(-t)*heaviside(t);
%     dy0 = 2;
% d)
%    p1 = [1,2,5];
%    p2 = [1,0];
%    y0 = 2;
%    dy0 = 0;
%    x_t = heaviside(t);
%    
%%      
[EC_3,Yf,Yn,Yf_inv,Yn_inv] = Resolucion_ED(p1,p2,y0,x_t,dy0);