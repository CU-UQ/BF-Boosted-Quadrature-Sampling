function S = FastRegVol(X, s, lambda)
%FastRegVol Fast reverse iterative volume sampling for linear regression
%
%S = FastRegVol(X, s, lambda) samples s rows from X, and returns the
%indices of those rows in S. This code is an implementation of Alg. 3 in 
%[Derezinski and Warmuth, JMLR 19, 2018]. Lambda is a regularization
%parameter. 

[n, d] = size(X);

Z = inv(X.'*X + lambda*eye(d));
S = 1:n;
while length(S) > max(s, 2*d)
    while true
        sample = S(randsample(length(S), 1));
        hi = 1 - X(sample, :)*Z*X(sample, :).';
        if binornd(1, hi)
            break
        end
    end
    S = setdiff(S, sample);
    Zx = Z*X(sample, :).';
    Z = Z + Zx*Zx.'/hi; 
end

if s < 2*d
    S = RegVol(X(S, :), s, lambda);
end

end

