function coefficients = levinson(corr)
% Input
% corr              autocorrelation vectors of size M + 1
% Output
% coefficients      linear prediction coefficients of size M
% | r(1)      r(2)        r(M-2)  r(M) | a(1)   |   | r(2)   |
% | r(2)      r(1)        r(M-2)  r(M) | a(2)   |   | r(3)   |
% |                                    |        | = |        |
% | r(M-1)    r(M-2)      r(1)    r(2) | a(M-1) |   | r(M)   |
% | r(M)      r(M-1)      r(2)    r(1) | a(M)   |   | r(M+1) |
M = length(corr) - 1;
E(1) = corr(1);
kappa(1) = corr(2);
a(1,1) = kappa(1);
E(2) = E(1)*(1-kappa(1)^2);

m = 2;
while (m <= M)
    m3 = m;
    corr3 = corr(m+1);
    as = a(1:1:m-1,m-1);
    coors = corr(m:-1:2);
    q(m) = corr(m+1) - sum(a(1:1:m-1,m-1)'.*corr(m:-1:2));
    kappa(m) = q(m)/E(m);
    a(m,m) = kappa(m);
    a(1:m-1,m) = a(1:m-1,m-1) - kappa(m)*a(m-1:-1:1,m-1);
    E(m+1) = E(m)*(1-kappa(m)^2);
    m = m + 1;
end
coefficients = a(:,M)';