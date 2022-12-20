function M = segmentation_cantidad_ciclos(y,m,t,t_m,cant_ciclos,N)
%m: es el vector de remuestreo que se usó para construir y
%t: es el eje temporal de la funcion que se remuestreó
%t_m: es el eje temporal de la función remuestrada
%cant_ciclos: cantidad de ciclos
%N: duracion de cada ciclo

[~,comienzo] = min(abs(t));

[~,A]  = min(abs(m-comienzo));

P = cant_ciclos;
M = zeros(P,N); %matriz que aloja a los periodos

for p = 1:P
    M(p,:) = interp1(0.5*(t_m-t_m(A))/pi,y,p-1+[0:1/(N):1-1/(N)],'spline');%acá decidimos que el muestreo de cada periodo tenga N muestras
end;
