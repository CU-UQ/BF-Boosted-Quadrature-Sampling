function [xn, wn, y_LF, y_HF, problem_dim] = data_loader(dataset)
%data_loader Auxiliary function used to load data
%
%[xn, wn, y_LF, y_HF] = data_loader(dataset) returns 1D quadrature points xn
%with corresponding weights wn, as well as low- and high-fidelity data in
%the vectors y_LF and y_HF, respectively. 

switch dataset
    case "2d_cavity_flow"
        load("cavity_flow_dataset/2d_cavity_flow/Uc_tp_grid_2d_size_16_N_40.mat", 'Uc')
        y_LF = mean(Uc, 1).';
        load("cavity_flow_dataset/2d_cavity_flow/Uc_tp_grid_2d_size_128_N_40.mat", 'Uc')
        y_HF = mean(Uc, 1).';
        N = 40;
        left_end_pt = -1;
        right_end_pt = 1;
        [xn, wn] = lgwt(N, left_end_pt, right_end_pt);
        problem_dim = 2;
        
    case "3d_cavity_flow"
        load("cavity_flow_dataset/3d_cavity_flow/Uc_tp_grid_3d_size_16_N_10.mat", 'Uc')
        y_LF = mean(Uc, 1).';
        load("cavity_flow_dataset/3d_cavity_flow/Uc_tp_grid_3d_size_128_N_10.mat", 'Uc')
        y_HF = mean(Uc, 1).';
        N = 10;
        left_end_pt = -1;
        right_end_pt = 1;
        [xn, wn] = lgwt(N, left_end_pt, right_end_pt);
        problem_dim = 3;
        
    case "6d_cavity_flow"
        load("cavity_flow_dataset/6d_cavity_flow/Uc_tp_grid_6d_size_16_N_5", 'Uc')
        y_LF = mean(Uc, 1).';
        load("cavity_flow_dataset/6d_cavity_flow/Uc_tp_grid_6d_size_128_N_5", 'Uc')
        y_HF = mean(Uc, 1).';
        N = 5;
        left_end_pt = -1;
        right_end_pt = 1;
        [xn, wn] = lgwt(N, left_end_pt, right_end_pt);
        problem_dim = 6;
        
    case "duffing"
        load("Datasets\duffing_oscillator\duffing_lf.mat", 'xn', 'wn', 'y_LF')
        load("Datasets\duffing_oscillator\duffing_hf.mat", 'y_HF')
        problem_dim = 2;
        
    case "exponential"
        N = 21;
        left_end_pt = -1;
        right_end_pt = 1;
        [y_HF, xn, wn] = generate_exp_data(N, left_end_pt, right_end_pt);
        y_LF = 1 + repmat(xn,21,1) + repelem(xn,21,1); % Linear approximation to exp(x+y) at (0,0)
        problem_dim = 2;
end

end
