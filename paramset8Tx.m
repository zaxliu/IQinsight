function [params] = paramset8Tx()
    %% File I\O
    params.FILE1 = 'bb_signal_8Tx_1frame.txt';         % Data source file
    %% OFDM Mod\Demod  
    params.nc = 8;                                     % Num of carriers
    params.nSlot = 20;                                 % Num of slots recorded
    params.nOfdmSym = 7;                               % Num of OFDM symbols per slot
    params.mode = [100 90 90 90 90 90 90];             % Num of sample points for CP          
    params.nFFT = 1280;                                % Num of sub-carriers (FFT points)
    params.idxsc = [641:1280,1:640].';                 % Sub-carrier index correction        
    %% Field & Version Mapping
    for x=1
    nF = params.nFFT; nT = params.nSlot*params.nOfdmSym;   
    params.nField = 88;                                % Num of global stores needed
    params.fieldMap = zeros(nF,params.nc,nT);          % Mapping of F,T,Ant to Field
    params.versionMap = zeros(nF,params.nc,nT);        % Mapping of F,T,Ant to Field Version
    params.occurField = zeros(params.nField,nT);
    params.occurVersion = cell(params.nField,nT);      % Which version occured in certain time
    params.gResetFlag = zeros(params.nField,nT);       % Whether to reset a field's global store in certain time
    a = 40; b = 41;                                    % Shift const. to translate position in 1200 to 1280
    for x = 1
    %%%% May seem overlap, but already taken care by the order    
    %%%%#1-40 Data x4 
    for iF = 1:40
        for iV = 1:4
            iiF = iT((1:60)+mod((iF-1),20)*60);
            iiA = iV + 4*floor((iF-1)/20);
            params.fieldMap(iiF,iiA,:) = iF;
            params.versionMap(iiF,iiA,:) = iV;
        end
    end
    params.gResetFlag(1:40,2:14:nT) = 1;
    %%%%#41-80 DMRS x4
    for iF = 41:80
        for iV = 1:4
            iiF2 = (1:60)+mod(iF-41,20)*60;
            iiF = iT(intersect(union(union(2:12:1200,7:12:1200),12:12:1200),iiF2));
            iiA = iV + 4*floor((iF-41)/20);
            iiT = union(6:7:nT,7:7:nT);
            params.fieldMap(iiF,iiA,iiT) = iF;
            params.versionMap(iiF,iiA,iiT) = iV;
            params.fieldMap(:,:,6);
        end
    end
    params.gResetFlag(41:80,6:14:nT) = 1;
    %%%%#81 Guard Band & Blank x1
    params.fieldMap([1:40,1242:end],:,:) = 81;                       % Guard band
    params.fieldMap([(565:600)+a,(601:636)+b],:,[1:14,71:84]) = 81;  % GSS
    params.versionMap(params.fieldMap==81) = 1;
    params.gResetFlag(81,1) = 1;
    %%%%#82 CI1 x4 
    params.fieldMap(iT(1:1200),1:4,1:14:end) = 82;
    for nV = 1:4
        params.versionMap(iT(1:1200),nV,1:14:end) = nV;
    end
    params.gResetFlag(82,1) = 1;
    %%%%#83 CI2 x4
    params.fieldMap(iT(1:1200),5:8,1:14:end) = 83;
    for nV = 1:4
        params.versionMap(iT(1:1200),nV+4,1:14:end) = nV;
    end
    params.gResetFlag(83,1) = 1;
    %%%%#84 PBCH x8    
    params.fieldMap([(565:600)+a,(601:636)+b],1:8,8:11) = 84;    
    for nV = 1:8
        params.versionMap([(565:600)+a,(601:636)+b],nV,8:11) = nV;
    end
    params.gResetFlag(84,8) = 1;
    %%%%#85 PSS x16 
    params.fieldMap([(570:600)+a,(601:631)+b],:,7) = 85;           % PSS1
    params.fieldMap([(570:600)+a,(601:631)+b],:,77) = 85;          % PSS2
    for nV = 1:8
        params.versionMap([(570:600)+a,(601:631)+b],nV,7) = nV;
    end
    for nV = 9:16
        params.versionMap([(570:600)+a,(601:631)+b],nV-8,77) = nV;
    end
    params.gResetFlag(85,[7,77]) = 1;
    %%%%#86 SSS1 x8
    params.fieldMap([(570:600)+a,(601:631)+b],:,6) = 86;
    params.fieldMap([(570:600)+a,(601:631)+b],:,76) = 86;
    for nV = 1:8
        params.versionMap([(570:600)+a,(601:631)+b],nV,[6,76]) = nV;
    end
    params.gResetFlag(86,[6,76]) = 1;

    %%%%#87 RS x8   
    params.fieldMap([(1:6:600)+a,(601:6:1200)+b],5:8,1:7:nT) = 87; % RS0 on ant 1
    params.fieldMap([(4:6:600)+a,(604:6:1200)+b],5:8,5:7:nT) = 87; % RS0 on ant 1
    params.fieldMap([(1:6:600)+a,(601:6:1200)+b],1:4,5:7:nT) = 87; % RS1 on ant 0
    params.fieldMap([(4:6:600)+a,(604:6:1200)+b],1:4,1:7:nT) = 87; % RS1 on ant 0
    for nV = 1:4
        params.versionMap([(1:6:600)+a,(601:6:1200)+b],nV+4,1:7:nT) = nV; % RS0 on ant 1
        params.versionMap([(4:6:600)+a,(604:6:1200)+b],nV+4,5:7:nT) = nV; % RS0 on ant 1
    end
    for nV = 5:8
        params.versionMap([(1:6:600)+a,(601:6:1200)+b],nV-4,5:7:nT) = nV; % RS1 on ant 0
        params.versionMap([(4:6:600)+a,(604:6:1200)+b],nV-4,1:7:nT) = nV; % RS1 on ant 0
    end
    params.gResetFlag(87,1) = 1;
%     params.gResetFlag(87,2) = 1;
    %%%%#81
    params.fieldMap([(1:6:600)+a,(601:6:1200)+b],1:4,1:7:nT) = 81;   % RS0 on ant 0
    params.fieldMap([(4:6:600)+a,(604:6:1200)+b],1:4,5:7:nT) = 81;   % RS0 on ant 0
    params.fieldMap([(1:6:600)+a,(601:6:1200)+b],5:8,5:7:nT) = 81;   % RS1 on ant 1
    params.fieldMap([(4:6:600)+a,(604:6:1200)+b],5:8,1:7:nT) = 81;   % RS1 on ant 1
    params.versionMap(params.fieldMap==81) = 1;
    %%%%#88 DC x1
    params.fieldMap(641,:,:) = 88;                         
    params.versionMap(641,:,:) = 1;
    params.gResetFlag(88,1) = 1;
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
    %% Compression Ratio Calculation
    params.IQBits = 30;
    params.constBits = 30;
    params.vqBits = zeros(88,1);
    params.vqBits(1:40) = 6;     % 64QAM
    params.vqBits(41:80) = 2;    % QPSK
    params.vqBits(81) = 0;       % =0
    params.vqBits(82:83) = 4;    % 9 points
    params.vqBits(84) = 2;       % +-1 + 0
    params.vqBits(85) = 4;       % 16 points
    params.vqBits(86) = 1;       % 2 points
    params.vqits(87) = 2;        % QPSK
    params.vqBits(88) = 0;       % Constant
    params.weightBits = 32;
    params.weightBitsIntSign = 4;
    %% EVM
    params.n_evmblock = 20;
    params.evmblocklen = 9600;
    %% HDL Test
    % Save a filed at a specifice symbol for HDL test
    params.testOFDMsym = 2;      % Symbol index of test version occurence 
    params.testField = 40;       % Field index of test version occurence
end
function [idx2] = iT(idx1)
    a = 40; b = 41; center = 600;
    idx2 = zeros(size(idx1));
    idx2(idx1<=center) = idx1(idx1<=center) + a;
    idx2(idx1>center) = idx1(idx1>center) + b;
end