function lev_scores = FastLevScoreEstimation(A, r1, r2)
%FastLevScoreEstimation Computes leverage score estimates via sketching
%
%lev_scores = FastLevScoreEstimation(A, r1, r2) takes an n-by-d
%tall-and-skinny (i.e., with n > d) input matrix A and returns a length-n
%vector lev_scores which contains estimates of the leverage scores of A.
%This function is an implementation of Algorithm 1 in [Drineas et al., JMLR
%13, 2012]. The parameters r1 and r2 control the sketching sizes used for
%the two sketches. We need r1 > n; see the paper by Drineas et al. for
%details.

[n, d] = size(A);

assert(r1 > d, 'r1 must be greater than d')
assert(n > d, 'A should be tall-and-skinny')

% Compute FJLT of A
sgn_vec = round(rand(n, 1)) * 2 - 1;
A_mixed = fft(sgn_vec .* A) / sqrt(n);
rand_samp = randsample(n, r1, true);
Pi1_A = sqrt(n / r1) * A_mixed(rand_samp, :);

% Compute SVD of SA and appropriate pseudo-inverse
[~, S, V] = svd(Pi1_A, 'econ');
R_pinv = V * pinv(S);

% Apply JLT to pseudo-inverse
Pi2 = sqrt(1/r2) * randn(d, r2);
R_pinv_Pi2 = R_pinv * Pi2;

% Compute estimate of (A * pinv(A)) and associated leverage score estimates
Omega = A * R_pinv_Pi2;
lev_scores = sum(abs(Omega).^2, 2);

end
