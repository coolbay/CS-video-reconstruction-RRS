%%   
% This program performs the CS image reconstruction method Phase1 - intra RRS  in 
% the paper "Video Compressive Sensing Reconstruction via Reweighted Residual
% Sparsity" published in T-CSVT. 


% Written by C Zhao, Jan. 2015.

clear
clc
cur = cd;
addpath(genpath(cur));

total_seq_num = 4;
total_frame_num = 17; % to be tuned
GOP_size = 8; % to be tuned

key_subrate = 0.7;
subrate = 0.2;  % uniform Sampling Rate

block_size = 32; % Block Size for BCS

% 
imgPSNR  = zeros(total_seq_num, total_frame_num);

for kk = 2   % for each sequence   modified 1:1
    switch kk
        case 1
            sequence_name = 'bus_cif.yuv';
        case 2
            sequence_name = 'foreman_cif.yuv';
        case 3
            sequence_name = 'akiyo_cif.yuv';
        case 4
            sequence_name = 'mobile_cif.yuv';
    end
    
    % Constructe Measurement Matrix (Gaussian Random)
    N = block_size * block_size;
    M_key = round(key_subrate * N);
    M = round(subrate * N);
    randn('seed',0);
    PhiN = orth(randn(N, N))';
    Phi_key = PhiN(1:M_key, :);
    Phi = PhiN(1:M, :);
    
    Opts = [];
    Opts.block_size = block_size;
    
    % encoder and initial-decoder process
    for i = 1 : total_frame_num    % for each frame of this sequence
        if  strcmp(sequence_name, 'susie')
            filename = ['..\Sequences\' sequence_name '_' num2str(i) '.pgm'];
        else
            filename = ['..\Sequences\' sequence_name '_' num2str(i) '.png'];
        end
        frame{i} = double(imread(filename));
        [row, col] = size(frame{i});
        Opts.filename = filename;
        Opts.row = row;
        Opts.col = col;
        
        if mod(i, GOP_size) == 1
            Opts.Phi = Phi_key;
            Opts.subrate = key_subrate;
        end
        if mod(i, GOP_size) ~= 1
            Opts.Phi = Phi;
            Opts.subrate = subrate;
        end
        y{i} = BCS_Encoder(frame{i}, Opts.Phi, Opts.block_size);
        tic
        [x_initial{i},psnr] = Intra_RRS(y{i}, frame{i}, Opts);
        toc
        imgPSNR(kk,i) = psnr;
    end
end
xlswrite('PSNR results_All sequences.xlsx',imgPSNR);
