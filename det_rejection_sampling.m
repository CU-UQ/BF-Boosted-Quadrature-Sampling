function S = det_rejection_sampling(X, q, k)
%det_rejection_sampling Implementation of leveraged volume sampling
%
%S = det_rejection_sampling(X, q, k) returns a subset S of k row indices of
%X, where q contains the leverage scores of X. This function does leveraged
%volume sampling via the determinantal rejection sampling algorithm in
%[Derezinski et al., arXiv:1802.06749, 2018].

[n, d] = size(X);
no_samp = max(k, 4*d^2);
q = q(:);  % Ensure it's a column vector
XTX_inv = inv(X.'*X);

accept = false;
while ~accept
    Pi = randsample(n, no_samp, true, q);
    X_Pi = 1./sqrt(q(Pi)) .* X(Pi, :);
    accept_prob = det((X_Pi.'*X_Pi)*XTX_inv/no_samp);
    if rand <= accept_prob
        accept = true;
    end
end

vol_S = FastRegVol(X_Pi, k, 0);

S = Pi(vol_S);

end
