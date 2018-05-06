%%   
% This program performs CS video reconstruction based on Phase 2 - inter RRS in 
% the paper "Video Compressive Sensing Reconstruction via Reweighted Residual
% Sparsity" published in T-CSVT. 


% Written by C Zhao, Jan. 2015.

clear
clc
cur = cd;
addpath(genpath(cur));

total_seq_num = 4;
total_frame_num = 17; 
GOP_size = 8; 

key_subrate = 0.7;
subrate = 0.2;  % uniform Sampling Rate

block_size = 32; % Block Size for BCS

% Constructe Measurement Matrix (Gaussian Random)
N = block_size * block_size;
M_key = round(key_subrate * N);
M = round(subrate * N);
randn('seed',0);
PhiN = orth(randn(N, N))';
Phi_key = PhiN(1:M_key, :);
Phi = PhiN(1:M, :);

% Record the PSNR results in an excel
imgPSNR  = zeros(total_seq_num, total_frame_num);

for j = 2
    switch j
        case 1
            sequence_name = 'bus_cif.yuv';
        case 2
            sequence_name = 'foreman_cif.yuv';
        case 3
            sequence_name = 'akiyo_cif.yuv';
        case 4
            sequence_name = 'mobile_cif.yuv';
    end
    
    disp('Initilization ...');
    
    % encoder and initial-decoder process
    for i = 1 : total_frame_num
        frame{i} = double(imread(['..\Sequences\' sequence_name '_' num2str(i) '.png']));             % read the orginal frame
        [row, col] = size(frame{i});
        
        % encode the orginal key frames and non-key frames
        if mod(i, GOP_size) == 1
            y{i} = BCS_Encoder(frame{i}, Phi_key, block_size);
        else
            y{i} = BCS_Encoder(frame{i}, Phi, block_size);
        end
        
        % read the initial recovery for the current frame
        Info_Dir = dir('..\Results_Phase1\');   % the directory that contains all the initial frames recovered from the first phase
        Dir_Num = length(Info_Dir);         % number of video frames in the directory
        Cmp_Name = strcat(sequence_name,'_',num2str(i));
        Com_Num = length(Cmp_Name);
        for kk=3:Dir_Num
            temp_Name = Info_Dir(kk).name;
            if strcmp(Cmp_Name,temp_Name(1:Com_Num))
                initial_name = Info_Dir(kk).name;
                break
            end
        end
        x_initial{i} = double(imread(['..\Results_Phase1\' initial_name]));
    end
    
    Opts = [];
    Opts.initial = x_initial;           % initial recovery for all the frames
    Opts.Phi_key = Phi_key;             % A
    Opts.Phi = Phi;
    Opts.row = row;
    Opts.col = col;
    Opts.max_iterations = 60;    
    Opts.Inloop = 300;
    Opts.mu = 2.5e-3;
    Opts.thr = 8;  
    Opts.org = frame;                   % all the original frames
    Opts.y = y;                         % the measurements  for all the frames
    Opts.block_size = block_size;
    Opts.frame_num = total_frame_num;
    
    
    disp('Beginning of inter RRS for video CS Recovery');
    % First process the key frames
    cur_no_array = [1,9,17,2,3,4,5,6,7,8,10,11,12,13,14,15,16];
    
    if ~exist('Convergence','dir') 
        mkdir('Convergence'); 
    end
    if ~exist('..\Results_Phase2\', 'dir')
        mkdir '..\Results_Phase2\'
    end
    % Refine the initial recovery using the temporal SLSM
    for i = 1 : total_frame_num %total_frame_num
        
        cur_no = cur_no_array(i);
        
        if cur_no == 1 || cur_no == 9|| cur_no == 17    
            final_recovery{cur_no} = Opts.initial{cur_no};
        else
            if cur_no < 9
                ref1 = 1;
                ref2 = 9;
            elseif cur_no < 17
                ref1 = 9;
                ref2 = 17;
            end          
            % the core function for inter RRS
            [frame_reconstructed, All_PSNR] = Inter_RRS(Opts,cur_no,ref1,ref2);
            final_recovery{cur_no} = frame_reconstructed;
        end
        
        psnr = PSNR(frame{cur_no}, final_recovery{cur_no});
        imgPSNR(j,cur_no) = psnr;
        
        
        Final_Name = strcat(sequence_name,'_',num2str(cur_no),'_PSNR_',num2str(csnr(frame{cur_no},final_recovery{cur_no},0,0)),'dB.png');
        imwrite(uint8(final_recovery{cur_no}),strcat('..\Results_Phase2\',Final_Name));
        
        Plot_flag = 1;
        if Plot_flag
            if  cur_no == 1 || cur_no == 9|| cur_no == 17
                disp('Because I am a key frame, I output nothing here!');
            else
                figure;imshow(uint8(final_recovery{cur_no}));title(['PSNR = ' num2str(psnr) ' dB']);
                disp(['PSNR = ' num2str(psnr) ' dB']);
                figure; plot(1:Opts.max_iterations, All_PSNR, 'LineWidth',2.0),
                title(strcat(sequence_name,'.00',num2str(cur_no),' subrate=',num2str(subrate),' Evolution of PSNR (dB)'));
                set(gca,'FontName','Times'),
                set(gca,'FontSize',14),
                xlabel('Iterative Numbers ');
                ylabel('Inter RRS');
                saveas(gcf,['Convergence\' Final_Name],'png');
                close all;
            end
        end       
    end    
    disp('End of Inter RRS');
    xlswrite('PSNR results_the inter RRS.xlsx',imgPSNR);
end
