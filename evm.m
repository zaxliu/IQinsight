function [EVM] = evm(bb_c,bb_c_rec,params)
    nc = params.nc;
    n_block = params.n_evmblock;
    BLOCK_LEN = params.evmblocklen;
    EVM = zeros(nc,n_block);
    for idx_block = 1:n_block
        %% IQ data load: pre-processing
        DATA_PRE = bb_c(:,((idx_block-1)*BLOCK_LEN+1):(idx_block*BLOCK_LEN));
        %% IQ data load:  post-processing
        DATA_POS = bb_c_rec(:,((idx_block-1)*BLOCK_LEN+1):(idx_block*BLOCK_LEN));

        %% Block EVM calculation
        error = DATA_POS - DATA_PRE;
        for idx =1:nc
            EVM(idx, idx_block) = sqrt((error(idx,:)*error(idx,:)')...
                /(DATA_PRE(idx,:)*DATA_PRE(idx,:)') )*100;
        end
    end

%% Output 
% disp('EVM(%) increment value for each carrier (row) and block (column):')
% disp(EVM)
disp('Averaged EVM(%) increment for each carrier:');
disp(mean(EVM.').');
disp('Max EVM(%) for each carrier:');
disp(max(EVM.').');
end