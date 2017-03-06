function [codebook,idxSQ] = coreSQ_li(r)
global params;
    [m,n] = size(r);
    r = reshape(r,m*n,1);
    %%% Partition    
    rr = real(r); mr = min(rr); Mr = max(rr); Ar = Mr - mr;
    ri = imag(r); mi = min(ri); Mi = max(ri); Ai = Mi - mi;
    % Split bits
    br = round(Ar^2/(Ar^2+Ai^2)*params.rrBits); bi = params.rrBits - br;
    if br == 0
        ir = ones(size(rr)).';
        codebookr = 0;
    else
        lenr = Ar/2^br;
        codebookr = linspace(mr,Mr,2^br);
        partitionr = (mr+lenr):lenr:(Mr-lenr/2);
        ir = quantiz(rr,partitionr) + 1;
    end
    if bi == 0
        ii = ones(size(ri)).';
        codebooki = 0;
    else
        leni = Ai/2^bi;
        codebooki = linspace(mi,Mi,2^bi);
        partitioni = (mi+leni):leni:(Mi-leni/2);
        ii = quantiz(ri,partitioni) + 1;
    end
    codebook{1} = codebookr;
    codebook{2} = codebooki;
    idxSQ{1} = ir;
    idxSQ{2} = ii;
end