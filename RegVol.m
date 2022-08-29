function S = RegVol(X, s, lambda)
%RegVol Reverse iterative volume sampling for linear regression
%
%S = RegVol(X, s, lambda) samples s rows from X, and returns the indices of
%those rows in S. This code is an implementation of Alg. 2 in [Derezinski
%and Warmuth, JMLR 19, 2018]. Lambda is a regularization parameter.

[n, d] = size(X);

Z = inv(X.'*X + lambda*eye(d));
ZX = Z*X.';
h = zeros(1,n);
for k = 1:n
    h(k) = 1 - X(k, :)*ZX(:, k);
end
S = 1:n;
while length(S) > s
    sample = S(randsample(length(S), 1, true, h(S)/sum(h(S))));
    S = setdiff(S, sample);
    v = (Z*X(sample, :).')/sqrt(h(sample));
    for k = 1:length(S)
        h(S(k)) = h(S(k)) - (X(S(k), :)*v)^2;
    end
    Z = Z + v*v.';
end

end

