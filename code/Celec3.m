function [elec]=Celec3()
global Vth boxes

elec(1,1)=rand()*200e-9;
elec(1,2)=rand()*100e-9;

    while(isbox(elec(1,1:2), boxes))
            elec(1,1:2) = [200e-9*rand 100e-9*rand];
    end

elec(1,3)=randn()*Vth/(2^0.5);
elec(1,4)=randn()*Vth/(2^0.5);

end