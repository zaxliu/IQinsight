function VerilogDataInit(data,fileDir)
%******************************************************
% data: input versions 
% fileDir: output file dir
% Versions:
%   YDB 12.4
%   LJC 12.4 Nested function dec2binPN
%******************************************************
[sx,sy]=size(data);
if(sx>800 || sy>8)
    disp('input data size error!')
else
    data_tmp = zeros(800,8);
    
    data_tmp(1:sx,1:sy)=data;
    for m = 1: 8
        I = round(real(data_tmp(:,m)));
        Q = round(imag(data_tmp(:,m)));
        if(length(fileDir)==0)
            filename = [fileDir,'version',num2str(m),'.coe'];
        else
            filename = [fileDir,'\version',num2str(m),'.coe'];
        end
        fid =fopen(filename,'w');
        for k=1:800
            tmp1 = dec2binPN(I(k),32);
            tmp2 = dec2binPN(Q(k),32);
            fprintf(fid,'%s',tmp1);
            fprintf(fid,'%s',tmp2);
            
            if(k<800)
                fprintf(fid,'\n');
            end
        end
        fclose(fid);
    end
end
end
function [numbin] = dec2binPN(numdec,N)
%判断输入数正负
if (numdec>= 0)
    %正数转二进制
    numbin1 = dec2bin(numdec,N);
else
    %负数转二进制
    numbin1 = dec2bin(abs(numdec),N);
    l1=length(numbin1);
    numbin4=0;
    for i=1:l1
        if (numbin1(l1-i+1)==num2str(1))%按位取反，用十进制表示
            numbin4=numbin4+0;
        else
            numbin4=numbin4+2^(i-1);
        end
    end
    %末位加1
    numbin4=numbin4+1;
    %把处理完的十进制数转成二进制，并输出在numbin
    numbin5=dec2bin(numbin4);
    numbin1=num2str(numbin5,N);
end
numbin=numbin1;
end
