function data = readVerilogData(fileDir,sx,sy)
%******************************************************
% fileDir: input verilog file dir
% data: verilog data :size:sx*sy
% sx,sy : actual data size
%******************************************************
if(sx>800 || sy>8)
    disp('input data size error')
else
    for m = 1:8
        if(length(fileDir)==0)
            filename = [fileDir,'Decom_',num2str(m),'.coe'];
        else
            filename = [fileDir,'\Decom_',num2str(m),'.coe'];
        end
        fid =fopen(filename,'r');
        tmp = fscanf(fid,'%d');
        fclose(fid);
        data_tmp(:,m)=tmp(1:2:end)+i*tmp(2:2:end);
    end
    data=data_tmp(1:sx,1:sy);
end

