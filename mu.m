function corr = mu(A,samp_list,b)
% This function generates optimality coefficient of the given sketching and
% QoI vector b.
% 
% Input:
% A: Matrix A of the sketched linear system SAx = Sb;
% samp_list: Sampling index of the sketching S;
% b: Vector b of the linear system SAx = Sb.

SQ = A(samp_list,:);
Q_perp = null(A');
y = Q_perp*Q_perp'*b;
Sy = y(samp_list);
x = pinv(SQ)*Sy;
corr = norm(x)/norm(y);
end