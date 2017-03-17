clc; clear; close all;

%% Variables de entrada

maxit = 10000;  %numero maximo de iteraciones a realizar

cr = 20;    %numero de caminos posibles
gn = 10;    %maxima cantidad de cambios de direccion

pc = 0.8;   %probabilidad de cruze
pm = 0.05;  %probabilidad de mutacion

cobst = 4;  %cantidad de obstaculos

limx = 500;     %ancho del area de busqueda
limy = 500;     %largo del area de busqueda
ladobs = 50;    %lado del obstaculo

%% Crear y delimitar los obstaculos

obx = (limx - ladobs/2)*rand(1, cobst); oby = (limy - ladobs/2)*rand(1, cobst); %posicion central de los obstaculos en x y y

for o = 1:length(obx)
    
    xlimit = [obx(o) - ladobs/2, obx(o) + ladobs/2];    %limites de los obstaculos en el plano para
    ylimit =  [oby(o) - ladobs/2, oby(o) + ladobs/2];   %x y y.
    
    xbox(:, o) = xlimit([1 1 2 2 1])';          %construccion de la funcion que describen  
    ybox(:, o) = ylimit([1 2 2 1 1])';          %los obstaculos en x y y.
    
end

%% Crar las principales variables de arranque

ang = pi*rand(gn + 1, cr);      %angulo de arranque 
vel = 90*rand(gn + 1, cr) + 10; %velocidades de arranque

%% Algoritmo genetico 

for l = 1:maxit
    clear c;    %se limpia en cada incio la funcion de costo
    
%%%%%%%%%%%%%%%%%%%%%%% Conversion de angulos y velocidades a rutas

    x(1, :) = ones(1, cr)*250;      %posiciones iniciales de x
    y(1, :) = ones(1, cr);              %posiciones iniciales de y
    
    for o = 1:cr        %construccion de las posiciones de acuerdo al angulo y velocidad
        for p = 1:gn + 1
            x(p+1, o) = x(p) + vel(p, o)*cos(ang(p, o));
            y(p+1, o) = y(p) + vel(p, o)*sin(ang(p, o));
        end
    end
    
%%%%%%%%%%%%%%%%%%%%%%% Obtencion del peso

r1 = sum(x < 0 | x > limx | y < 0 | y > limy) == 0;     %verifica cuales caminos sales del area de búsqueda
r2 = find(r1==0);   %asigna a r2 las rutas que se pierden
r3 = find(r1);          %asigna a r3 las rutas que se mantienen en el area

for p = r3
    for o = 1:length(obx)
        [xi, ~] = polyxpoly(x(:,p) , y(:, p) , xbox(:,o), ybox(:,o));   %busca si existen cruces entre el camino y los obstaculos
        xin(p, o) = length(xi)/2;   %debido a que siempre hay una entrada y salida, solo se cuentan los choques
    end                                    %que se representan por cada dos choques
end

aux1 = sum(xin, 2);     %contabiliza la cantidad de obstaculos con los que interfiere
c = 7 - (5/(1 + max(aux1)))*aux1;   %dependiendo del maximo de obstaculos, evalua cada ruta dado su rendimiento con una calificacion maxima de 7.                                           

c(r3) = c(r3) + 3 - (3./(max(limy - y(end, r3)) + 1))*(limy - y(end,r3))';   %los ultimos 3 puntos se asignan dependiendo de la distancia a la que se encuentren del objetivo.

c(r2) = 0;  %asigna la calificacion de 0 a aquellos vectores que se perdieron en la ruta.

if sum(c > 9) > 0
    break;      %si se obtiene un puntuaje perfecto (no choques y termina en el objetivo) se finaliza la busqueda.
end

%%%%%%%%%%%%%%%%%%%%%%% Seleccion natural

NPx = ang(:, abs(c) >= abs(0.99*mean(c)));     %obtiene solo aquellos angulos mayores al promedio
NPy = vel(:, abs(c) >= abs(0.99*mean(c)));      %obtiene solo aquellas velocidades mayores al promedio

%%%%%%%%%%%%%%%%%%%%%%% Reproduccion

NPx1 = zeros(size(ang));    %prepara las variables que contendran los hijos 
NPy1 = zeros(size(vel));

ind2 = randi(size(NPx, 2), 2, cr);  %realiza las mezclas de padres y madres que se reproduciran,

for o = 1:2:cr
        
        n1 = ind2(find(c(ind2(:, o)) == max(c(ind2(:, o))), 1), o);     %busca quien de los dos padres es el mejor
        n2 = ind2(find(c(ind2(:, o + 1)) == max(c(ind2(:, o + 1))), 1), o + 1); %busca quien de las dos madres es mejor

        g = randi(gn-1);         % punto de corte notas
        
        NPx1(:, o) = [NPx(1:g, n1); NPx(g + 1:end, n2)];    %se realiza el cruce para la variable "angulos"
        NPx1(:, o + 1) = [NPx(1:g, n2); NPx(g + 1:end, n1)];
        
        NPy1(:, o) = [NPy(1:g, n1); NPy(g + 1:end, n2)];    %se realiza el cruce para la variable "velocidad"
        NPy1(:, o + 1) = [NPy(1:g, n2); NPy(g + 1:end, n1)];
        
end

%%%%%%%%%%%%%%%%%%%%%%% Mutacion

for o = 1:cr            %se evalua para cada cromosoma o ruta
        if rand <= pm %la probabilidad de si sufrira mutacion o no
            for u = 1:randi(gn);  %con el numero random de genes a mutar
                NPx1(randi(gn), o) = pi*rand;   %para las variables de angulo
                NPy1(randi(gn), o) = 90*rand + 10; %y velocidad de manera aleatoria
            end
        end
end

%%%%%%%%%%%%%%%%%%%%%%% Actualizacion de estados

ang = NPx1;     %se actualizan los angulos 
vel = NPy1;     %se actualizan las velocidades

for o = 1:cr        %se reconstruyen los nuevos caminos consiguiendo los puntos nuevos x y y.
    for p = 1:gn + 1
        x(p+1, o) = x(p) + vel(p, o)*cos(ang(p, o));
        y(p+1, o) = y(p) + vel(p, o)*sin(ang(p, o));
    end
end

%%%%%%%%%%%%%%%%%%%%%%% Gráfica de resultado

clf     %se limpia el espacio 
hold on     

for o = 1:length(obx)       %se plotean cada uno de los obstaculos
    mapshow(xbox(:, o) ,ybox(:, o) ,'DisplayType', 'polygon', 'LineStyle', 'none', 'FaceColor', 'r')
end

for j=1:size(x, 2)          %se plotean cada uno de las rutas
    plot(x(:,j),y(:,j))
end

axis([0 500 0 500])     %se delimita el area de busqueda
hold off
pause(eps)     %se permite observar el resultado gráficamente

end