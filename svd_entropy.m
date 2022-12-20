function e = svd_entropy(M)

[~,s,~] = svd(M);
S = diag(s);
S = S/sum(S);
e = -sum(S.*log(S));