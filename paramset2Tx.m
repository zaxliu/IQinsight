function [params] = paramset2Tx()
    %% File I\O
    params.FILE1 = 'bb_signal_2Tx_1frame.txt';         % Data source file
    %% OFDM Mod\Demod  
    params.nc = 2;                                     % Num of carriers
    params.nSlot = 20;                                 % Num of slots recorded
    params.nOfdmSym = 7;                               % Num of OFDM symbols per slot
    params.mode = [160 144 144 144 144 144 144];             % Num of sample points for CP          
    params.nFFT = 2048;                                % Num of sub-carriers (FFT points)
    params.idxsc = [1025:2048,1:1024].';                 % Sub-carrier index correction        
    %% Field & Version Mapping
    for x=1
    nF = params.nFFT; nT = params.nSlot*params.nOfdmSym;   
    params.nField = 10;                                % Num of global stores needed
    params.fieldMap = zeros(nF,params.nc,nT);          % Mapping of F,T,Ant to Field
    params.versionMap = zeros(nF,params.nc,nT);        % Mapping of F,T,Ant to Field Version
    params.occurField = zeros(params.nField,nT);
    params.occurVersion = cell(params.nField,nT);      % Which version occured in certain time
    params.gResetFlag = zeros(params.nField,nT);       % Whether to reset a field's global store in certain time
    for x = 1
    %%%% May seem overlap, but already taken care by the order    
    %%%%#1-2 Data x1 
    for iF = 1:2
        iiF = 1:nF;
        iiA = iF;
        params.fieldMap(iiF,iiA,:) = iF;
        params.versionMap(iiF,iiA,:) = 1;
    end
    params.gResetFlag(1,2) = 1;
    params.gResetFlag(2,2) = 1;
    %%%%#3 Guard Band & idl & Blank x1
    params.fieldMap([1:424,1626:end],:,:) = 3;                       % Guard band + idl
    params.fieldMap(iT(565:636),:,[1:14,71:84]) = 3;  % GSS
    params.versionMap(params.fieldMap==3) = 1;
    params.gResetFlag(3,1) = 1;
    %%%%#4 CI1 x1 
    params.fieldMap(iT(1:1200),1,1:14:end) = 4;
    params.versionMap(params.fieldMap==4) = 1;
    params.gResetFlag(4,1) = 1;
    %%%%#5 CI2 x1
    params.fieldMap(iT(1:1200),2,1:14:end) = 5;
    params.versionMap(params.fieldMap==5) = 1;
    params.gResetFlag(5,1) = 1;
    %%%%#6 PBCH x2    
    params.fieldMap(iT(565:636),1:2,8:11) = 6;    
    params.versionMap(params.fieldMap==6) = 1;
    params.gResetFlag(6,8) = 1;
    %%%%#7 PSS x4 
    params.fieldMap(iT(570:631),:,7) = 7;           % PSS1
    params.fieldMap(iT(570:631),:,77) = 7;          % PSS2
    for nV = 1:2
        params.versionMap(iT(570:631),nV,7) = nV;
    end
    for nV = 3:4
        params.versionMap(iT(570:631),nV-2,77) = nV;
    end
    params.gResetFlag(7,[7,77]) = 1;
    %%%%#8 SSS1 x4
    params.fieldMap(iT(570:631),:,6) = 8;           % SSS1
    params.fieldMap(iT(570:631),:,76) = 8;          % SSS2
    for nV = 1:2
        params.versionMap(iT(570:631),nV,[6,76]) = nV;
    end
    params.gResetFlag(8,[6,76]) = 1;
    %%%%#9 RS x2   
    params.fieldMap(iT(1:6:1200),2,1:7:nT) = 9; % RS0 on ant 1
    params.fieldMap(iT(4:6:1200),2,5:7:nT) = 9; % RS0 on ant 1
    params.fieldMap(iT(1:6:1200),1,5:7:nT) = 9; % RS1 on ant 0
    params.fieldMap(iT(4:6:1200),1,1:7:nT) = 9; % RS1 on ant 0
    params.versionMap(iT(1:6:1200),2,1:7:nT) = 1; % RS0 on ant 1
    params.versionMap(iT(4:6:1200),2,5:7:nT) = 1; % RS0 on ant 1
    params.versionMap(iT(1:6:1200),1,5:7:nT) = 2; % RS1 on ant 0
    params.versionMap(iT(4:6:1200),1,1:7:nT) = 2; % RS1 on ant 0
    params.gResetFlag(9,1) = 1;
    params.gResetFlag(9,2) = 1;
    %%%%#3
    params.fieldMap(iT(1:6:1200),1,1:7:nT) = 3;   % RS0 on ant 0
    params.fieldMap(iT(4:6:1200),1,5:7:nT) = 3;   % RS0 on ant 0
    params.fieldMap(iT(1:6:1200),2,5:7:nT) = 3;   % RS1 on ant 1
    params.fieldMap(iT(4:6:1200),2,1:7:nT) = 3;   % RS1 on ant 1
    params.versionMap(params.fieldMap==3) = 1;
    %%%%#10 DC x1
    params.fieldMap(1025,:,:) = 10;                         
    params.versionMap(1025,:,:) = 1;
    params.gResetFlag(10,1) = 1;
    end
    %%%% occurField
    for idxT = 1:nT
        for idxF = 1:params.nField
            params.occurField(idxF,idxT) = sum(sum(params.fieldMap(:,:,idxT)==idxF))~=0;
        end
    end
    %%%% occurVersion
    for idxT = 1:nT
        field = params.fieldMap(:,:,idxT);
        version = params.versionMap(:,:,idxT);
        for idxF = 1:params.nField
            params.occurVersion{idxF,idxT} = ...
                unique(version(field==idxF));
        end
    end
    end
   %% VQ SQ
    params.VQth = 2000;
    params.rrBits = 8;
   %% Compression Ratio
    params.IQBits = 30;
    params.constBits = 30;       % Bitwidth for constellation codebook
    params.vqBits = zeros(88,1);
    params.vqBits(1:2) = 2;      % QPSK
    params.vqBits(3) = 0;        % =0
    params.vqBits(4:5) = 4;      % 9 points
    params.vqBits(6) = 2;        % +-1,0
    params.vqBits(7) = 4;        % 16 points
    params.vqBits(8) = 1;        % 2 points
    params.vqits(9) = 2;         % QPSK
    params.vqBits(10) = 0;       % Constant
    params.weightBits = 32;
    params.weightBitsIntSign = 4;
    %% EVM
    params.n_evmblock = params.nSlot;
    params.evmblocklen = sum(params.mode)+params.nOfdmSym*params.nFFT;
   %% HDL Test
    % Save a filed at a specifice symbol for HDL test
    params.testOFDMsym = 1;      % Symbol index of test version occurence 
    params.testField = 4;        % Field index of test version occurence
end
function [idx2] = iT(idx1)
    a = 424; b = 425; center = 600;
    idx2 = zeros(size(idx1));
    idx2(idx1<=center) = idx1(idx1<=center) + a;
    idx2(idx1>center) = idx1(idx1>center) + b;
end