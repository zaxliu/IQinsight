function [idxVQ,valVQ] = coreVQ(refVersion,gIdx)
    global gConst gNumConst params;
    num = length(refVersion);
    idxVQ = zeros(num,1);
    valVQ = zeros(num,1);
    for n = 1:num
        if gNumConst(gIdx)==0  %% can do before for-loop, but need context info
            gConst{gIdx}(1,1) = refVersion(n);
            gNumConst(gIdx) = 1;
            idxVQ(n) = 1;
            valVQ(n) = refVersion(n);
        else
            [v,i] = min(abs(gConst{gIdx}(1:gNumConst(gIdx))-refVersion(n)));
            if v>params.VQth
                gNumConst(gIdx) = gNumConst(gIdx)+1;
                gConst{gIdx}(gNumConst(gIdx),1) = refVersion(n);
                idxVQ(n) = gNumConst(gIdx);
                valVQ(n) = refVersion(n);
            else
                idxVQ(n) = i;
                valVQ(n) = gConst{gIdx}(i);
            end
        end
    end
end
