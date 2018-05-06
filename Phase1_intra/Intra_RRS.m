%
%
function [reconstructed_image,psnr] = Intra_RRS(y, original_image, Opts)

Phi = Opts.Phi;
subrate = Opts.subrate;
row  = Opts.row;
col = Opts.col;
block_size = Opts.block_size;
filename = Opts.filename;

Thr_Value = 8;

% Obtain Initilization by MH
disp('Initilization ...');
[x_MH x_DWT] = MH_BCS_SPL_Decoder(y, Phi, subrate, row, col);

if ~isfield(Opts,'initial')
    Opts.initial = double(x_MH);
end

Opts.org = original_image;

if ~isfield(Opts,'IterNum')
    Opts.IterNum = 40; %120
end

if ~isfield(Opts,'thr')
    Opts.thr = Thr_Value;
end

if ~isfield(Opts,'mu')
    Opts.mu = 2.5e-3;
end

if ~isfield(Opts,'Inloop')
    Opts.Inloop = 300; % 300
end

fprintf('Initial PSNR = %0.2f\n',csnr(Opts.org,Opts.initial,0,0));
% Invoke Proposed LDCT Alogorithm for Block-based CS Recovery
disp('Beginning of LDCT Algorithm for CS Recovery');
[reconstructed_image All_PSNR]= SBI_Iter(y, Opts);

psnr = PSNR(original_image, reconstructed_image);

if ~exist('Convergence\', 'dir')
    mkdir 'Convergence\'
end

if ~exist('..\Results_Phase1\', 'dir')  
    mkdir '..\Results_Phase1\'
end

save(['Convergence\' filename '_rate_' num2str(subrate) '_results.mat'], 'All_PSNR');

Final_Name_LDCT = strcat(filename,'_rate_',num2str(subrate),'_thr_',num2str(Thr_Value),'_LDCT_SBI_Iter','_PSNR_',num2str(csnr(original_image,reconstructed_image,0,0)),'dB.tif');
imwrite(uint8(reconstructed_image),strcat('Results_Phase1\',Final_Name_LDCT));

Plot_flag = 1;
if Plot_flag
    figure(100);imshow(uint8(reconstructed_image));title(['LDCT PSNR = ' num2str(psnr) ' dB']);
    disp(['LDCT PSNR = ' num2str(psnr) ' dB']);
    figure; plot(1:Opts.IterNum,All_PSNR, 'LineWidth',2.0),
    title(strcat(filename,' subrate=',num2str(subrate),' Evolution of PSNR (dB)'));
    set(gca,'FontName','Times'),
    set(gca,'FontSize',14),
    xlabel('Iterative Numbers ');
    ylabel('PSNR');
    saveas(gcf,['Convergence\' Final_Name_LDCT],'png');
    close all;
end
disp('End of LDCT');

