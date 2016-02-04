function r = autocorr(input)
% Calculate autocorrelation of input signal
% Input
%   (vector)    input   windows of N samples
% Output
%   (vector)    r       correlation values
N = length(input)
s = input(1:N);
r = zeros(1, N);
for i = 1:N
    st(1:N) = zeros(1, N);
    st(1:N+1-i) = input(i:N);
    r(i) = sum(s.*st);
end