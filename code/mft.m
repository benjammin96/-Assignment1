function [MFT]=mft(JXY)
global NumP 

for i=1:NumP
    N=nnz(JXY(i,:));
    MFT(i,1)=sum(JXY(i,:)/N);
end


end