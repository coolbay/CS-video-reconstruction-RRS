
function [x_reconstructed, All_PSNR] = BCS_LDCT_Decoder_Bframe(Opts, cur,ref1,ref2)

row = Opts.row;
col = Opts.col;
Phi = Opts.Phi;
mu = Opts.mu;
Inloop = Opts.Inloop;
max_iterations = Opts.max_iterations;
x_initial = Opts.initial{cur};
block_size = Opts.block_size;
y = Opts.y{cur};
x_org = Opts.org{cur};

Para = [];
Para.v = Opts.thr;

x_ref1 = Opts.initial{ref1};
x_ref2 = Opts.initial{ref2};

x = im2col(x_initial, [block_size block_size], 'distinct');
b = zeros(size(x));

All_PSNR = zeros(1,max_iterations);

ATy = Phi'*y;
ATA = Phi'*Phi;
IM = eye(size(ATA));

for i = 1:max_iterations
   
    x_hat = x;
    % first solve Eq. (16), with z_j being x_hat and b_j being ...b   (b_0 being 0)
    r = col2im(x_hat - b, [block_size block_size],[row col], 'distinct');  % x_hat: z_j
    
    % Construct groups; for each group, recover each patch; use
    % multi-hypothesis to recover the entire image x_j+1
    x_bar = LDCT_Solver_Bframe(r, x_ref1, x_ref2, Para);
    
    x_bar = im2col(x_bar, [block_size block_size], 'distinct');
    
    u = x_bar; % u:x_j+1
    
    % then solve Eq. (15), with x_j+1 being u, b_j being... b
    % and z_j being x_hat
    for kk = 1:Inloop       
        d = ATA*x_hat - ATy + mu*(x_hat - u - b);
        dTd = d'*d;
        G = d'*(ATA + mu*IM)*d;
        Step_Matrix = abs(dTd./G);
        Step_length = diag(diag(Step_Matrix));
        x = x_hat - d*Step_length;
        x_hat = x;        
    end
    % x_hat, x: z_j+1
    
    b = b - (x - u); % x: z_j+1; u:x_j+1 (其作用只在于计算一个b)
    
    x_img = col2im(x, [block_size block_size],[row col], 'distinct');
    
    Cur_PSNR = csnr(x_img,x_org,0,0);
    All_PSNR(i) = Cur_PSNR;
    fprintf('IterNum = %d, PSNR = %0.2f\n',i,Cur_PSNR);
end
x_reconstructed = x_img;



function [ImgRec] = LDCT_Solver_Bframe(x, x_ref1, x_ref2, Para)

if ~isfield(Para,'PatchSize')
    Para.PatchSize = 8;
end

if ~isfield(Para,'Profile')
    Para.Profile = 'fast';
end

if strcmp(Para.Profile, 'normal')
    SlidingDis   =  2;
elseif strcmp(Para.Profile, 'fast')
    SlidingDis   =  4;
end

if ~isfield(Para,'Factor')
    Para.Factor = 2.8*6;
end

if ~isfield(Para,'ArrayNoCur')
    Para.ArrayNoCur = 1;
end

if ~isfield(Para,'ArrayNoRef1')
    Para.ArrayNoRef1 = 5;
end

if ~isfield(Para,'ArrayNoRef2')
    Para.ArrayNoRef2 = 5;
end

if ~isfield(Para,'SearchWin')
    Para.SearchWin = 20;
end

Threshold = Para.Factor*Para.v;


[Hight, Width]   =    size(x);
SearchWin        =    Para.SearchWin;
PatchSize        =    Para.PatchSize;
PatchSize2       =    PatchSize*PatchSize;
ArrayNoCur       =    Para.ArrayNoCur;
ArrayNoRef1      =    Para.ArrayNoRef1;
ArrayNoRef2      =    Para.ArrayNoRef2;
AllArrayNo = ArrayNoCur + ArrayNoRef1 + ArrayNoRef2;
Dict = DCT2D_Matrix(PatchSize)';

N     =  Hight-PatchSize+1;
M     =  Width-PatchSize+1;
L     =  N*M;

Row     =  [1:SlidingDis:N];
Row     =  [Row Row(end)+1:N];
Col     =  [1:SlidingDis:M];
Col     =  [Col Col(end)+1:M];

PatchSet      =  zeros(PatchSize2, L, 'single');
PatchSetRef1  =  zeros(PatchSize2, L, 'single');
PatchSetRef2  =  zeros(PatchSize2, L, 'single');

Count  =  0;
for i  = 1:PatchSize
    for j  = 1:PatchSize
        Count   =  Count+1;
        Patch   =  x(i:Hight-PatchSize+i,j:Width-PatchSize+j);
        Patch1  =  x_ref1(i:Hight-PatchSize+i,j:Width-PatchSize+j);
        Patch2  =  x_ref2(i:Hight-PatchSize+i,j:Width-PatchSize+j);
        Patch   =  Patch(:);
        Patch1  =  Patch1(:);
        Patch2  =  Patch2(:);
        PatchSet(Count,:)     =  Patch';
        PatchSetRef1(Count,:) =  Patch1';
        PatchSetRef2(Count,:) =  Patch2';
    end
end

PatchSetT     =   PatchSet';
PatchSetRef1T =   PatchSetRef1';
PatchSetRef2T =   PatchSetRef2';

I        =   (1:L);
I        =   reshape(I, N, M);
NN       =   length(Row);
MM       =   length(Col);

ImgTemp     =  zeros(Hight, Width);
ImgWeight   =  zeros(Hight, Width);
PatchArray  =  zeros(PatchSize, PatchSize, ArrayNoCur);

%tic;

for  i  =  1 : NN
    for  j  =  1 : MM
        
        CurRow      =   Row(i);
        CurCol      =   Col(j);
        Off         =   (CurCol-1)*N + CurRow;
        
        v = PatchSetT(Off, :);
        CurPatchIndx    =  PatchSearchValue(PatchSetT, CurRow, CurCol, v, ArrayNoCur, SearchWin, I);
        CurPatchIndx(1) = Off;
        Ref1PatchIndex  =  PatchSearchValue(PatchSetRef1T, CurRow, CurCol, v, ArrayNoRef1, SearchWin, I);
        Ref2PatchIndex  =  PatchSearchValue(PatchSetRef2T, CurRow, CurCol, v, ArrayNoRef2, SearchWin, I);
        
        CurArray  = PatchSet(:, CurPatchIndx);
        Ref1Array = PatchSetRef1(:, Ref1PatchIndex);
        Ref2Array = PatchSetRef2(:, Ref2PatchIndex);
        
        hp = 80;
        Array = [CurArray Ref1Array Ref2Array];
        %[SS UU VV] = svd(Array);
        %Dict = SS';
        DifArray = Array - repmat(v',1, AllArrayNo);
        Dis = sum(DifArray.^2)/(PatchSize2);
        wei         =  exp( -Dis./hp );
        Weight         =  wei./(sum(wei)+eps);
        
        MX = sum(repmat(Weight,PatchSize2,1).* Array,2);
        
        alphaArray = Dict * Array;
        beta = Dict * MX;
        betaArray = repmat(beta,1,AllArrayNo);
        RX = Dict*(Array - repmat(MX,1,AllArrayNo));
        SO     =   mean(RX.^2, 2);
        SO     =   max(0, SO-Para.v);
        tau   =   (0.56*sqrt(2)*Para.v^2)./(sqrt(SO) + eps);
        tauArray = repmat(tau,1,AllArrayNo);
        
        CurArray    =   Dict'*(soft( alphaArray-betaArray, tauArray ) + betaArray);
        
        
        for k = 1:ArrayNoCur
            PatchArray(:,:,k) = reshape(CurArray(:,k),PatchSize,PatchSize);.
        end
        
        for k = 1:length(CurPatchIndx)
            RowIndx  =  ComputeRowNo((CurPatchIndx(k)), N);
            ColIndx  =  ComputeColNo((CurPatchIndx(k)), N);
            ImgTemp(RowIndx:RowIndx+7, ColIndx:ColIndx+7)    =   ImgTemp(RowIndx:RowIndx+7, ColIndx:ColIndx+7) + PatchArray(:,:,k)';
            ImgWeight(RowIndx:RowIndx+7, ColIndx:ColIndx+7)  =   ImgWeight(RowIndx:RowIndx+7, ColIndx:ColIndx+7) + 1;
        end
        
    end
end


ImgRec = ImgTemp./(ImgWeight+eps);



function  INDX  =  PatchSearchValue(X, Row, Col, v, Nv, S, I)
[N, M]    =   size(I);
Dim2      =   size(X,2);

rmin    =   max( Row-S, 1 );
rmax    =   min( Row+S, N );
cmin    =   max( Col-S, 1 );
cmax    =   min( Col+S, M );

idx     =   I(rmin:rmax, cmin:cmax);
idx     =   idx(:);
B       =   X(idx, :);

dis     =   (B(:,1) - v(1)).^2;

for k = 2:Dim2
    dis   =  dis + (B(:,k) - v(k)).^2;
end
dis   =  dis./Dim2;

[~, ind]   =  sort(dis);

INDX        =  idx( ind(1:Nv) );













