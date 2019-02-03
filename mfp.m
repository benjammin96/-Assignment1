function [MFP]=mfp(MFPX, MFPY)
global  MFPL MFP N MaxIt NumP 

for j=1:MaxIt
    for i =1:NumP
        MFPL(i,j)=sqrt(MFPX(i,j)^2 +MFPY(i,j)^2);
    end
end

for i=1:NumP
    N=nnz(MFPL(i,:));
    MFP(i,1)=sum(MFPL(i,:)/N);
end


end