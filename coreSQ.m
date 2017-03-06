function [codebook,idxSQ] = coreSQ(r)
global params;
    [m,n] = size(r);
    r = reshape(r,m*n,1);
    rr = real(r); ri = imag(r);
    %%% Load Partition and Codebook 
    b = params.rrBits/2;
    if b==0
        [codebook,idxSQ] = coreSQ_li(r);
    else
        codebookName = ['lloydCodebook' num2str(b) 'bit.mat'];
        load(codebookName);
        % Split bits
        ir = quantiz(rr,partitionr) + 1;
        ii = quantiz(ri,partitioni) + 1;
        codebook{1} = codebookr;
        codebook{2} = codebooki;
        idxSQ{1} = ir;
        idxSQ{2} = ii;
    end
end