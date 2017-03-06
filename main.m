%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Script for BB Compression
% Version:
%   LJC 9.29
%   LJC 11.4
%   LJC 11.9
%   LJC 12.2 Change defination of compression ratio
%   YDB 12.4 X and Y label for plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
global gConst gNumConst gWeight gRefVer params;
%% Initializtion
%%% Parameters
% params = paramset2Tx;                               % Set parameters
params = paramset8Tx;                               % Set parameters
fid1 = fopen(params.FILE1); frewind(fid1);          % Data File
nOfdmSym = params.nOfdmSym;                         % Num of OFDM symbols per slot
nOfdmSymTotal = params.nSlot*nOfdmSym;              % Total num of OFDM symbols                         
mode = params.mode;                                 % Num of sample points for CP & OFDM symbols
nc = params.nc;                                     % Num of carriers
nFFT = params.nFFT;                                 % Num of sub-carriers (FFT points)

%%% Global Variables
gConst = cell(params.nField,1);                     % (global) constellation points
gNumConst = zeros(params.nField,1);                 % (global) Num of constellation points found
gWeight = cell(params.nField,1);                    % (global) Versions weight
gRefVer = cell(params.nField,1);                    % (global) Reference version (version 1)

%%% Test Related Variables
fff = zeros(nFFT,nc,140);                          % F-domain original signal
% ttt = zeros(1280,8,140);                          % T-domain original signal
fff_rec = zeros(nFFT,nc,140);                      % F-domain recovered signal
% ttt_rec = zeros(1280,8,140);                      % T-domain recovered signal
bbi = 0;
bbri = 0;
bb = zeros(nc,(sum(mode)+length(mode)*nFFT)*params.nSlot); 
bb_rec = zeros(nc,(sum(mode)+length(mode)*nFFT)*params.nSlot);
% time_comp = zeros(1,nOfdmSymTotal);
% time_decomp = zeros(1,nOfdmSymTotal);
% rrr = zeros(nFFT,nc,140);

bVQSQ = zeros(nOfdmSymTotal,1);
bWeight = zeros(nOfdmSymTotal,1);
bConstCodebook = zeros(nOfdmSymTotal,1);

if exist('versions.mat','file')
    delete versions.mat
end

%% Compression
resultC = cell(2,nOfdmSymTotal);
for idx = 1:nOfdmSymTotal
%     clc; disp(idx);
%     tic;
    %%% Load data and convert to freq domain (OFDM demodulation)    
    nCP = mode(mod(idx-1,7)+1);                         % Position in slot
    tt = fscanf(fid1,'%X\n',2*nc*(nCP+nFFT));           % Read data from file
    tt(tt>32767) = tt(tt>32767) - 65536;                % Convert to signed int
    tt = tt(1:2:end) + 1j*tt(2:2:end);                  % Convert to complex    
    tt = reshape(tt,nc,(nFFT+nCP));                     % Realign by nFFT x nc
    bb(:,(bbi+1):(bbi+nCP+nFFT)) = tt;
    bbi = bbi + nCP + nFFT;
    tt = tt(:,(nCP+1):end);                             % Discard CP
    tt = tt.';
    ff = fft(tt,nFFT,1);                                % Get freq domain modulation symbol 
    ff = ff(params.idxsc,:);
    fff(:,:,idx) = ff;

    %%% Get weight of versions with respect to reference version and perform vector quantization
    resultC_temp1 = cell(3,length(params.nField));   % Compression result in current slice
    resultC_temp2 = cell(2,1);
    rr = zeros(size(ff));                            % Residual error
    gResetFlag = params.gResetFlag(:,idx);           % Reset flag for fields
    % Iterate through fields to perform VQ, This part can ran in parrallel
    for iF = 1:params.nField;                        
        if params.occurField(iF,idx)==0
            continue;                                % Bypass unappeared fields
        end
        %%%% Get data of different versions
        occurVersion = params.occurVersion{iF,idx};
        lenVer = sum(sum(...
            (params.fieldMap(:,:,idx)==iF)&...
            (params.versionMap(:,:,idx)==occurVersion(1))...
            ));
        versions = zeros(lenVer,length(occurVersion));
        for iV = 1:length(occurVersion)
            versions(:,iV) = ...
                ff((params.fieldMap(:,:,idx)==iF)...
                &(params.versionMap(:,:,idx)==occurVersion(iV)));
        end
        %%%% Save user specified field for HDL test
        if idx==params.testOFDMsym && iF==params.testField
            save('versions.mat','versions');
        end
        %%%% Vector quantization and compression        
        if gResetFlag(iF)==1                % Reset global store when new field appear
            gConst{iF} = [];
            gNumConst(iF) = 0;
        end
        if sum(occurVersion==1)~=0          % When new reference version show up...
            oldNumConst = gNumConst(iF);
            [idxVQ,valVQ] = coreVQ(versions(:,1),iF);   % Perform VQ (refresh constellation)
            [idxVQ_c] = coreLosslC(idxVQ);              % Perform lossless compression
            gRefVer{iF} = valVQ;                        % Store ref ver. in case other version show up in following epochs
            diffNumConst = gNumConst(iF)-oldNumConst;
            if diffNumConst~=0                          % Only transmit new constellation points
                resultC_temp1{2,iF} = ...
                    gConst{iF}((gNumConst(iF)-diffNumConst+1):gNumConst(iF),1);       
            end
            resultC_temp1{3,iF} = idxVQ_c; 
        end
        
        %%%%  Get weight of appeared version                                                                   !!!bit width        
        if gResetFlag(iF)==1                        % Reset global store when new version appear
            weight = coreWeight(versions,iF);       % Get version weight (Assume stored in float point fommat)
            gWeight{iF}(occurVersion) = weight;     % Store weight for computation that follows    
            resultC_temp1{1,iF} = weight;            
        else                                        
            weight = gWeight{iF}(occurVersion);     % Use stored weight
        end
               
        %%%%  Get residual error
        rec = floor(gRefVer{iF} * weight);                 % ??? Bitwidth ???
        for iV = 1:length(occurVersion)            
            rr(...
                (params.fieldMap(:,:,idx)==iF)...
                &(params.versionMap(:,:,idx)==occurVersion(iV))...
                ) = versions(:,iV) - rec(:,iV);
        end
        
        %%%% Record compression ratio related information
        bWeight(idx) = bWeight(idx) + params.weightBits * (length(occurVersion)-1) * gResetFlag(iF);
        bConstCodebook(idx) = bConstCodebook(idx) + params.constBits * diffNumConst;
        bVQSQ(idx) = bVQSQ(idx) + ...
            params.rrBits * length(occurVersion)*lenVer + ...
            params.vqBits(iF)* lenVer;
    end

    %%% Scaler quantization for residual Error
%   rrr(:,:,idx) = rr;
    [codebook,idxSQ] = coreSQ(rr);
    resultC_temp2{1} = codebook;
    resultC_temp2{2} = idxSQ;
    
    %%% Transmit compression result
    resultC{1,idx} = resultC_temp1;  
    resultC{2,idx} = resultC_temp2;
    
    %%% Time consumption
%    time_comp(idx) = toc;
end
%% Decompression
%%% Global Variables
gConst = cell(params.nField,1);                     % (global) constellation points
gNumConst = zeros(params.nField,1);                 % (global) Num of constellation points found
gWeight = cell(params.nField,1);                    % (global) Versions weight
gRefVer = cell(params.nField,1);                    % (global) Reference version (version 1)

for idx = 1:nOfdmSymTotal
%     clc; disp(idx);
%     tic;
    %%% Receive compression result
    resultC_temp1 = resultC{1,idx};
    resultC_temp2 = resultC{2,idx};
    
    %%% Recover residual error
    codebook = resultC_temp2{1};
    idxSQ = resultC_temp2{2};
    [rr] = codebook{1}(idxSQ{1})+1j*codebook{2}(idxSQ{2});
    rr = reshape(rr.',nFFT,nc);
    
    %%% Recover VQ and multiply weight    
    gResetFlag = params.gResetFlag(:,idx);   
    ff = zeros(nFFT,nc);
    for iF = 1:params.nField
        if params.occurField(iF,idx)==0
            continue;
        end        
        occurVersion = params.occurVersion{iF,idx};
        
        %%%%  VQ and lossless compression recovery
        if gResetFlag(iF)==1                % Reset global store when new field appear
            gConst{iF} = [];
            gNumConst(iF) = 0;
        end
        if sum(occurVersion==1)~=0
            if ~isempty(resultC_temp1{2,iF})
                gConst{iF} = [gConst{iF};resultC_temp1{2,iF}];
            end
            idxVQ_c = resultC_temp1{3,iF};
            [idxVQ] = coreLosslC_rec(idxVQ_c);
            gRefVer{iF} = gConst{iF}(idxVQ);
        end
        
        %%%%  Get weight
        if gResetFlag(iF)==1            
            weight = resultC_temp1{1,iF};            
            gWeight{iF}(occurVersion) = weight;            
        else                                        
            weight = gWeight{iF}(occurVersion);
        end
        
        %%%%  Get recovered data, save for HDL test and map to F-T grid
        lenVer = sum(sum(...
            (params.fieldMap(:,:,idx)==iF)&...
            (params.versionMap(:,:,idx)==occurVersion(1))...
            ));
        versions_rec = zeros(lenVer,length(occurVersion));
        for iV = 1:length(occurVersion)
            versions_rec(:,iV) =  rr((params.fieldMap(:,:,idx)==iF)&(params.versionMap(:,:,idx)==occurVersion(iV)))...
                + floor(gRefVer{iF} * weight(iV));
            ff((params.fieldMap(:,:,idx)==iF)&(params.versionMap(:,:,idx)==occurVersion(iV))) =...
                versions_rec(:,iV);               
        end
        if idx==params.testOFDMsym && iF==params.testField
            save('versions.mat','versions_rec','-append');
        end
    end
    
    %%% convert to time domain and save data
%     fff_rec(:,:,idx) = ff;
    ff(params.idxsc,:) = ff;
    tt = ifft(ff,nFFT);                             % Get freq domain modulation symbol
    tt = tt.';
    nCP = mode(mod(idx-1,7)+1);                     
    bb_rec(:,(bbri+1):(bbri+nCP+nFFT)) = [tt(:,(nFFT-nCP+1):nFFT),tt];
    bbri = bbri + nCP + nFFT;
    
    %%% Time consumption
%   time_decomp(idx) = toc;
end

%% Performance Evaluation
%%% EVM
[EVM] = evm(bb,bb_rec,params);

%%% Time Consumption
% plot(1:140,time_comp,1:140,time_decomp);

%%% Compression Ratio
bSym = params.IQBits * nc * (repmat(mode',params.nSlot,1) + nFFT);
bFrame = params.IQBits * nc * (sum(mode)+length(mode)*nFFT) * params.nSlot;
ratio = (bVQSQ+bWeight+bConstCodebook)./bSym;
ratioVQSQ = bVQSQ./bSym;
ratioWeight = bWeight./bSym;
ratioConst = bConstCodebook./bSym;
microAvg = mean(ratio);
macroAvg = sum(bVQSQ+bWeight+bConstCodebook)/bFrame;
varCoef = var(ratio)/microAvg^2;
 
plot(ratio,'b');hold on;
plot(ones(length(ratio),1)*microAvg,'b.');
plot(ratioVQSQ,'r');
plot(ones(length(ratioVQSQ),1)*mean(ratioVQSQ),'r.');
legend('Sym-wise Ratio w/  codebook','Micro Average w/ codebook','Sym-wise Ratio w/o codebook','Micro Average w/o codebook');
xlabel('OFDM Symbol Index');ylabel('Compression Ratio')
grid on
%
disp('W/  Codebook:');
disp(['Macro Average Compression Ratio:   ' num2str(macroAvg)]);
disp(['Micro Average Compression Ratio:   ' num2str(microAvg)]);
disp(['Variance Coefficient:   ' num2str(varCoef)]);
%
disp('W/O codebook:');
disp(['Macro Average Compression Ratio:   ' num2str(sum(bVQSQ)/bFrame)]);
disp(['Micro Average Compression Ratio:   ' num2str(mean(ratioVQSQ))]);
disp(['Variance Coefficient:   ' num2str(var(ratioVQSQ)/mean(ratioVQSQ).^2)]);
%% Post Processing
fclose(fid1);