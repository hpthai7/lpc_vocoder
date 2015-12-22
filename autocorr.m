function r = autocorr(input)
% input windows of N items
% vector of correlation, r(k), k = 1,K
N = length(input);
s = input(1:N);
r = zeros(1, N);
for i = 1:N
    st(1:N) = zeros(1, N);
    st(1:N+1-i) = input(i:N);
%     s
%     st
    r(i) = sum(s.*st);
end