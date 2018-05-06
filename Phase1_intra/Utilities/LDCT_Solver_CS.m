function [ImgRec] = LDCT_Solver_CS(ImgInput, Opts)

if ~isfield(Opts,'PatchSize')
    Opts.PatchSize = 8;
end

if ~isfield(Opts,'Profile')
    Opts.Profile = 'fast';
end

if strcmp(Opts.Profile, 'normal')
    SlidingDis   =  2;
elseif strcmp(Opts.Profile, 'fast')
    SlidingDis   =  4;
end

if ~isfield(Opts,'Factor')
    Opts.Factor = 2.8*6;
end

if ~isfield(Opts,'ArrayNo')
    Opts.ArrayNo = 10;
end

if ~isfield(Opts,'SearchWin')
    Opts.SearchWin = 20;
end


[Hight Width]   =   size(ImgInput);
SearchWin = Opts.SearchWin;
PatchSize    =    Opts.PatchSize;
PatchSize2    =   PatchSize*PatchSize;
ArrayNo   =   Opts.ArrayNo;
v = Opts.v;
Threshold = Opts.Factor*v;
Dict = DCT2D_Matrix(PatchSize)';
N     =  Hight-PatchSize+1;
M     =  Width-PatchSize+1;
L     =  N*M;

Row     =  [1:SlidingDis:N];
Row     =  [Row Row(end)+1:N];
Col     =  [1:SlidingDis:M];
Col    =  [Col Col(end)+1:M];


PatchSet     =  zeros(PatchSize2, L, 'single');

Count     =  0;
for i  = 1:PatchSize
    for j  = 1:PatchSize
        Count    =  Count+1;

        Patch  =  ImgInput(i:Hight-PatchSize+i,j:Width-PatchSize+j);
        Patch  =  Patch(:);

        PatchSet(Count,:) =  Patch';
    end
end

PatchSetT  =   PatchSet';

I        =   (1:L);
I        =   reshape(I, N, M);
NN       =   length(Row);
MM       =   length(Col);

ImgTemp     =  zeros(Hight, Width);
ImgWeight   =  zeros(Hight, Width);
IndcMatrix  =  zeros(NN, MM, ArrayNo);
PatchArray  =  zeros(PatchSize, PatchSize, ArrayNo);


for  i  =  1 : NN
    for  j  =  1 : MM
        
        CurRow      =   Row(i);
        CurCol      =   Col(j);
        Off      =   (CurCol-1)*N + CurRow;

        [CurPatchIndx Weight]  =  PatchMatching(PatchSetT, CurRow, CurCol, Off, ArrayNo, SearchWin, I);
        
        CurPatchIndx(1) = Off;
                      
        CurArray = PatchSet(:, CurPatchIndx); 
        
        MX = sum(repmat(Weight',PatchSize2,1).* CurArray,2); 
        
        alphaArray = Dict * CurArray;   
        beta = Dict * MX;   
        betaArray = repmat(beta,1,ArrayNo);
        RX = Dict*(CurArray - repmat(MX,1,ArrayNo)); 
        SO     =   mean(RX.^2, 2);        
        SO     =   max(0, SO-v);
        tau   =   (0.56*sqrt(2)*v^2)./(sqrt(SO) + eps);  
        tauArray = repmat(tau,1,ArrayNo);
        
        CurArray    =   Dict'*(soft( alphaArray-betaArray, tauArray ) + betaArray);
              
        for k = 1:ArrayNo
            PatchArray(:,:,k) = reshape(CurArray(:,k),PatchSize,PatchSize);
        end
               
        for k = 1:length(CurPatchIndx)
            RowIndx  =  ComputeRowNo((CurPatchIndx(k)), N);
            ColIndx  =  ComputeColNo((CurPatchIndx(k)), N);
            ImgTemp(RowIndx:RowIndx+7, ColIndx:ColIndx+PatchSize-1)    =   ImgTemp(RowIndx:RowIndx+PatchSize-1, ColIndx:ColIndx+PatchSize-1) + PatchArray(:,:,k)';
            ImgWeight(RowIndx:RowIndx+7, ColIndx:ColIndx+PatchSize-1)  =   ImgWeight(RowIndx:RowIndx+PatchSize-1, ColIndx:ColIndx+PatchSize-1) + 1;
        end
        
    end
end

ImgRec = ImgTemp./(ImgWeight+eps);

return;



