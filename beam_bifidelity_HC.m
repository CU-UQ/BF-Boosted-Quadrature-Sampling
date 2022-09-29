% Load bifidelity data
N = 10;
D = 4;
left_end_pt = -1;
right_end_pt = 1;
[xn, wn] = lgwt(N, left_end_pt, right_end_pt);
load('Y.mat');
y_LF = Y_LF;
y_HF = Y_HF;
J = 2;

% Compute matrix A^(n) which is the same for all n
An = zeros(N, J+1);
for i = 1:N
    An(i, :) = my_legendre_1d(J,xn(i))';
end
An = sqrt(wn).*An;

% Compute leverage scores of A
[Un, ~] = qr(An, 0);
Un_sq = Un.^2;

% Compute big matrix A for tensor product space
A_ = An;
Wn = sqrt(wn);
for d = 2:D
    A_ = kron(A_, An);
    Wn = kron(Wn,sqrt(wn));
end

% Truncate A into hyperbolic cross space
sub = sub_tp_idx_set(D,J,'tensor product');
s = prod(sub+1,1);
l = 1:(1+J)^D;
ind = l(s <= J+1);
A = A_(:,ind);
[row,col] = size(A);
[U,~] = qr(A,0);
p = sum(U.^2,2)./col;

% Adjust the weight for y_LF and y_HF
yL = Wn.*y_LF;
yH = Wn.*y_HF;

% Exact solution
my_tic = tic;
x_sol = A \ yH;
err_exact = norm(A*x_sol - yH)/norm(yH);
toc_exact = toc(my_tic);


% Sampling setup
no_samp = 18;
no_trials = 100;
t_alev = zeros(no_trials,1);
resamp = 10;
corr_lev = zeros(no_trials,2);
corr_lev_vol = zeros(no_trials,2);
corr_unif = zeros(no_trials,2);
err_all_boosted = zeros(no_trials,resamp);
err_lev_score = zeros(no_trials,1);
err_lev_score_boosted = zeros(no_trials,1);
err_unif = zeros(no_trials,1);
err_unif_boosted = zeros(no_trials,1);
err_lev_vol = zeros(no_trials,1);
err_lev_vol_boosted = zeros(no_trials,1);
y = yH;

% QR sampled solution
qr_time = ceil(no_samp/col);
samp_list = ones(col*qr_time,1);
A_temp = A;
for i = 1:qr_time
    [~, ~, perm_vec] = qr(A_temp.', 0);
    samp = perm_vec(1:col);
    samp_list((i-1)*col+1:i*col) = samp';
    A_temp(samp,:) = [];
end
samp = samp_list(1:no_samp);
x_sol = A(samp,:) \ yH(samp);
err_QR = norm(A*x_sol - yH)/norm(yH);

 % Leverage score sampling 
for tr = 1:no_trials
    [x_sol, ~,samp_list] = hyperbolic_cross_sampling(N,D,J,An,Un_sq,no_samp,y,'alev');
    corr_lev(tr,1) = mu(A,samp_list,yL);
    corr_lev(tr,2) = mu(A,samp_list,yH);
    err_lev_score(tr) = norm(y-A*x_sol)/norm(y);

 % Boosted leverage score sampling
    best_err = Inf;
    for smp = 1:resamp
        [x_sol,~,samp_list, scaling] = hyperbolic_cross_sampling(N,D,J,An,Un_sq,no_samp,yL,'alev');
        current_err = norm(A*x_sol - yL)/norm(yL);
        err_all_boosted(tr, smp) = current_err;
        if current_err < best_err
            best_samp = samp_list;
            best_scaling = scaling;
            best_err = current_err;
        end
    end
    SA = A(best_samp, :);
    Sy = yH(best_samp);
    x_sol = SA \ Sy;
    err_lev_score_boosted(tr) = norm(A*x_sol - yH)/norm(yH);
 
 % Uniform sampling
    [x_sol, ~,samp_list] = hyperbolic_cross_sampling(N,D,J,An,Un_sq,no_samp,y,'unif');
    corr_unif(tr,1) = mu(A,samp_list,yL);
    corr_unif(tr,2) = mu(A,samp_list,yH);
    err_unif(tr) = norm(y-A*x_sol)/norm(y);

 % Boosted leverage score sampling
    best_err = Inf;
    for smp = 1:resamp
        [x_sol,~,samp_list, scaling] = hyperbolic_cross_sampling(N,D,J,An,Un_sq,no_samp,yL,'unif');
        current_err = norm(A*x_sol - yL)/norm(yL);
        err_all_boosted(tr, smp) = current_err;
        if current_err < best_err
            best_samp = samp_list;
            best_scaling = scaling;
            best_err = current_err;
        end
    end
    SA = A(best_samp, :);
    Sy = yH(best_samp);
    x_sol = SA \ Sy;
    err_unif_boosted(tr) = norm(A*x_sol - yH)/norm(yH);
    
 % Leveraged volume sampling 
    samp = det_rejection_sampling(A, p, no_samp);
    scaling = 1./sqrt(p(samp));
    SA = scaling.*A(samp, :);
    Sy = scaling.*yH(samp);
    x_sol = SA \ Sy; 
    corr_lev_vol(tr,1) = mu(A,samp,yL);
    corr_lev_vol(tr,2) = mu(A,samp,yH);
    err_lev_vol(tr) = norm(A*x_sol - yH)/norm(yH);

 % Boosted leveraged volume sampling
        best_err = Inf;
        for smp = 1:resamp
            samp = det_rejection_sampling(A, p, no_samp);
            scaling = 1./sqrt(p(samp));
            SA = scaling.*A(samp, :);
            Sy = scaling.*yL(samp);
            x_sol = SA \ Sy; 
            current_err = norm(A*x_sol - yL)/norm(yL);
            if current_err < best_err
                best_samp = samp;
                best_err = current_err;
            end
        end
        SA = A(best_samp, :);
        Sy = yH(best_samp);
        x_sol = SA \ Sy;   
        err_lev_vol_boosted(tr) = norm(A*x_sol - yH)/norm(yH);
    tr
end

% save computation results
save('B_no_18_HC.mat','err_exact','err_QR','err_lev_score','err_lev_score_boosted','err_unif','err_unif_boosted','err_lev_vol','err_lev_vol_boosted');
save('B_Cor_no_18_HC.mat','corr_lev','corr_unif','corr_lev_vol');
 