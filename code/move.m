function[electrons] = move(electrons)
global dt Tmn eu Vth
    a= electrons(:,2)+electrons(:,4)*dt >100E-9; % logical indexing for if the electron goes above the top limit
    b =electrons(:,2)+electrons(:,4)*dt<0; %logical indexing for if the electron goes below the bottom limit
    electrons(:,4)=electrons(:,4) - 2*electrons(:,4).*a -2*electrons(:,4).*b  ; % Electron speed is equal to electron speed unless a logical indexing flags as a 1, then it reverses (relflects it back)
    electrons(:,1)=mod(electrons(:,1)+electrons(:,3)*dt,200E-9); % electron x- position is equal to the modulus of electron x- position with 200E-9 so anything below 200E-9 will just be regular but if greater than 200E-9 it will start back at the start
    electrons(:,2)=electrons(:,2)+electrons(:,4)*dt; % electron y co-ordinate is equal to electron y co-ordinate + speed
end