function [U] = my_istct_mod(F,f,g)
if nargin<3
   g = 0.3;
end

f = f-1+eps;

[N, ~] = size(F);

if iscolumn(f)
   f = f';
end

E = exp(-2*pi*1i*repmat(f,[N 1]).*repmat(1./f',[1 N])); 

aux = abs(F).^g;
aux = aux - repmat(mean(aux),[N,1]);

U = E * aux;
end