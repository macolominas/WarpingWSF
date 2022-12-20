function M = segmentation(y,t_m,cant_ciclos,N)


P = cant_ciclos;
M = zeros(P,N); 

for p = 1:P
    M(p,:) = interp1(0.5*(t_m)/pi,y,p-1+[0:1/(N):1-1/(N)],'spline');
end;
