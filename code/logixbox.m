function [alogic,ax1] = logixbox(electrons)
 %checks if in box
 global NumP  boxes dt  ax1
   
 alogic=zeros(NumP,1);
    
    for i=1:size(boxes,1)
        for j= 1:NumP  
            A=(electrons(j,1)+electrons(j,3)*dt*2);
            B=(electrons(j,1)+electrons(j,3)*dt*2);
            C=(electrons(j,2)+electrons(j,4)*dt*2);
            D=(electrons(j,2) + electrons(j,4)*dt*2);
            if( A> boxes(i,1) && B < boxes(i,2) && C > boxes(i,3) && D < boxes(i,4))
                alogic(j,1)=1;
                Delta1= abs(electrons(j,1)-boxes(i,1));
                Delta2= abs(electrons(j,1)- boxes(i,2));
                ax1(j,1)=Delta1 < 1E-9 || Delta2 <1E-9;
            end
        end 
    end
end