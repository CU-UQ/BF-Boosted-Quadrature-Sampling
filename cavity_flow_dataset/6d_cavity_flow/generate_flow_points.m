% Generate inputs that can be put into the flow code

N = 5;
D = 6;
left_end_pt = -1;
right_end_pt = 1;

xi = zeros(N^D, 52);

% Compute Gauss-Legendre quadrature points and weights
[xn, wn] = lgwt(N, left_end_pt, right_end_pt);

for d = 1:D
    xi(:,d) = repmat(repelem(xn, N^(d-1), 1), N^(6-d), 1);
end

save('xi_6d_cavity_flow', 'xi')
