clear;
close all;
%You must make sure that you do not click off of figure(2) before the movie
%is over or else it will plot over the histogram
global Vth dt j Tmn eu Pscat NumP MFPX XSUM YSUM MFPY MFPL MFP MaxIt JSUM JXY
m0 =9.11E-31; % Electron rest mass
mn=0.6*m0; % Effective Electron mass
kB=1.3806E-23; % Boltzmann Constant
T=300; % Temperature
Vth=(kB*T/mn)^0.5; % Thermal Velocity
Tmn=0.2E-12; % Time between collisions
Dmn=0.2E-12*Vth;% Mean Free Path
NumP = 1000; % Number of particles
MaxIt =100; % Maximum Iterations
ylimit=100E-9;
xlimit=200E-9;
dt= ylimit/(Vth*100); % time step 
color=hsv(5);
eu=2.71828;
temp=zeros(1,MaxIt);
Pscat = zeros(NumP,1);
Pscat(:,1)= 1 - (eu^(-1*dt/Tmn));

XSUM =zeros(NumP,1);
MFPX=zeros(NumP,MaxIt);
YSUM =zeros(NumP,1);
MFPY=zeros(NumP,MaxIt);
MFPL=zeros(NumP,MaxIt);
MFP=zeros(NumP,1);

JSUM=zeros(NumP,1);
JXY=zeros(NumP,MaxIt);




for j=1:NumP % creates electrons
      electrons(j,:)=Celec2();
      POS1(j,:)=POS(electrons);% Tracks initial position
      POS2(j,:)=POS(electrons);%Tracks initial position
end
    
MB = sqrt(electrons(:,3).^2 + electrons(:,4).^2);
MBAv = sum(MB)/length(MB);
figure(1)
histogram(MB);
title('Histogram of Electron Speeds');
xlabel('Speed (m/s)');
ylabel('Number of particles');
text(1.3E5,90,sprintf('Average Speed%d',MBAv)) % displays the temperature of the system on the plot
pause(1)



%Construct the frame.
figure(2)
xlim([0 200E-9]);
ylim([0 100E-9]);
hold on    
for j=1:MaxIt
    electrons = move2(electrons); % moves the electrons
    temp(j) = (sum(electrons(:,3).^2) + sum(electrons(:,4).^2))*mn/(kB*NumP); % Temperature for the system, do I need to divide by the number of electrons? 
    for i=1:NumP
          if abs(POS1(i,1) - electrons(i,1))>100E-9; % If the change between old and new is very large we jumped across the screen
            POS1(i,1)=electrons(i,1); %just say the old position is also the new position to avoid plotting across the screen
          end
          if(i<5) % only plots the first 5 electrons
             plot([POS1(i,1) electrons(i,1)],[POS1(i,2) electrons(i,2)],'Color',color(i,:)) %plots the electron line between old and new.
             pause(0.0001)
          end
    end
    POS1(:,1)=electrons(:,1);
    POS1(:,2)=electrons(:,2);

end

figure(2)
%makes the plot titles for moving electrons
title('electron trajectory')
xlabel('x (m)')
ylabel('y (m)')

%makes the temperature plot 
tx = [1 MaxIt];
figure(3)
plot(tx,temp(tx));
title('Temperature over Time')
xlabel('Time (10^-14 s)')
ylabel('Temperature (K)')
text(130E-9,100E-9,sprintf('Temperature =%d',temp(j))) % displays the temperature of the system on the plot

MFP =mfp(MFPX,MFPY);
MFT=mft(JXY);

figure(4)
plot(MFP,'*');
hold on
line([0 1000],[Dmn,Dmn],'Color','red','LineStyle','-');
title('Mean Free Path of Each Electron')
xlabel ('Electrons')
ylabel('Individual Electron Mean Free Path')
xlim([0 NumP]);
ylim([0 200E-9]);
text(200,150E-9,sprintf('Red Line = Theoretical Mean Free Path %d',Dmn)) % displays the temperature of the system on the plot
text(200,130E-9,sprintf('Blue Markers are average measured Mean Free Paths (m)'))

figure(5)
plot(MFT,'*');
hold on
line([0 1000],[Tmn,Tmn],'Color','red','LineStyle','-');
title('Mean Free Time of Each Electron')
xlabel ('Electrons')
ylabel('Individual Electron Mean Free Time (s)')
xlim([0 NumP]);
ylim([0 20E-13]);
text(200,1.8E-12,sprintf('Red Line = Theoretical Mean Free time %d',Tmn)) % displays the temperature of the system on the plot
text(200,1.6E-12,sprintf('Blue Markers are average measured Mean Free Times'))

