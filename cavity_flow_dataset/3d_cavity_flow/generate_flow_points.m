% Generate inputs that can be put into the flow code

N = 10;
D = 3;
left_end_pt = -1;
right_end_pt = 1;

xi = zeros(N^D, 52);

% Compute Gauss-Legendre quadrature points and weights
[xn, wn] = lgwt(N, left_end_pt, right_end_pt);

xi(:,1) = repmat(xn, N^2, 1);
xi(:,2) = repmat(repelem(xn, N, 1), N, 1);
xi(:,3) = repelem(xn, N^2, 1);

save('xi_3d_cavity_flow', 'xi')