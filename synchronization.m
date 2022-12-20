function FINAL = synchronization(M)

[K,N] = size(M);
FINAL = zeros(K,N);
E = zeros(2*K,2*K);

for i = 1:K

    E(2*(i-1)+1:2*i,2*(i-1)+1:2*i) = [1 0; 0 1];
    for j = i+1:K
        d = distancias(M(i,:),M(j,:));
        [~,theta] = min(d);
        theta = theta - 1;

        R = [cos(2*pi*theta/N) sin(2*pi*theta/N); -sin(2*pi*theta/N) cos(2*pi*theta/N)];

        E(2*(i-1)+1:2*i,2*(j-1)+1:2*j) = R;
        E(2*(j-1)+1:2*j,2*(i-1)+1:2*i) = R';
    end;
end;


[U,~,~] = svd(E);

Q = U(:,1:2);
R = cell(1,K);

for i = 1:K
    [PHI,~,PSI] = svd(Q(2*(i-1)+1:2*i,:));
    R{i} = PHI*PSI';
end;


for i = 1:K

    d = sign(R{i}(1,2))*round(N*(acos(R{i}(1,1))/(2*pi)));
    FINAL(i,:) = circshift(M(i,:),-d);
end;
d = sign(R{1}(1,2))*round(N*(acos(R{1}(1,1))/(2*pi)));
FINAL = circshift(FINAL,d,2);



