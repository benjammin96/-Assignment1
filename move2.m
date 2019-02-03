function[electrons] = move(electrons)
global dt Pscat NumP Vth MFPX XSUM j YSUM MFPY JSUM JXY
   
    
    a= electrons(:,2)+electrons(:,4)*dt >100E-9; % logical indexing for if the electron goes above the top limit
    b =electrons(:,2)+electrons(:,4)*dt<0; %logical indexing for if the electron goes below the bottom limit
    sp=rand(NumP,1);
    c = Pscat(:,1) >sp;
    
    JXY(:,j)=JSUM.*c*dt;
    JSUM=(JSUM+1).*~c;
    
    XSUM=XSUM + electrons(:,3).*dt.*~c;
    MFPX(:,j)=abs(XSUM).*c;

    YSUM = YSUM +electrons(:,4).*dt.*~c;
    MFPY(:,j)=abs(YSUM).*c;
    
   
    electrons(:,4)=electrons(:,4).*~c - 2*electrons(:,4).*a -2*electrons(:,4).*b + randn()*(Vth/sqrt(2)).*~a.*~b.*c ; % Logical indexing for top and bottoms and also scattering condition
    electrons(:,3)=electrons(:,3).*~c + randn()*(Vth/sqrt(2)).*~a.*~b.*c; % Scattering condution for x 
    electrons(:,1)=mod(electrons(:,1)+electrons(:,3)*dt,200E-9); % electron x- position is equal to the modulus of electron x- position with 200E-9 so anything below 200E-9 will just be regular but if greater than 200E-9 it will start back at the start
    electrons(:,2)=electrons(:,2)+electrons(:,4)*dt; % electron y co-ordinate is equal to electron y co-ordinate + speed
end