function weight = coreWeight(versions,iF)
global gRefVer params;
    weight = versions(1,:)./gRefVer{iF}(1);
    A = params.weightBits/2;    % Total bits for real or imag part
    a = A-params.weightBitsIntSign;  % Bits for fractional part, (1bit for sign, rest bits for integer part)
    weight = floor(weight*2^a)/2^a;
end