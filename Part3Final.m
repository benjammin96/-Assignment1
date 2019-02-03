clear;
close all;

%You must make sure that you do not click off of figure(2) before the movie
%is over or else it will plot over the histogram
global Vth dt j Tmn eu Pscat NumP  MaxIt boxes specular ax1 ax diffusive

specular=1;
diffusive=0;

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
eu=2.71828;
temp=zeros(1,MaxIt);
Pscat = zeros(NumP,1);
Pscat(:,1)= 1 - (eu^(-1*dt/Tmn));
ax1=zeros(NumP,1);
ax=zeros(NumP,1);

figure(2)
hold on
boxes = 1e-9.*[80 120 0 40; 80 120 60 100];

for j=1:size(boxes,1)
           plot([boxes(j, 1) boxes(j, 1) boxes(j, 2) boxes(j, 2) boxes(j, 1)],...
               [boxes(j, 3) boxes(j, 4) boxes(j, 4) boxes(j, 3) boxes(j, 3)], 'k-');
        end


for j=1:NumP % creates electrons
      electrons(j,:)=Celec3();
      POS1(j,:)=POS(electrons);% Tracks initial position
      POS2(j,:)=POS(electrons);%Tracks initial position
end
    


%Construct the frame.
figure(2)
xlim([0 200E-9]);
ylim([0 100E-9]);
hold on    
NumPP=20; % # of electrons to plot
color=hsv(NumPP);

for j=1:MaxIt
    electrons = move3(electrons); % moves the electrons
    temp(j) = (sum(electrons(:,3).^2) + sum(electrons(:,4).^2))*mn/(kB*NumP); % Temperature for the system, do I need to divide by the number of electrons? 
    for i=1:NumP
          if abs(POS1(i,1) - electrons(i,1))>100E-9; % If the change between old and new is very large we jumped across the screen
            POS1(i,1)=electrons(i,1); %just say the old position is also the new position to avoid plotting across the screen
          end
          if(i<NumPP) % only plots the first NumPP electrons
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


EDensity(:,1) = electrons(:,1);
EDensity(:,2) = electrons(:,2);
figure(3)
hist3(EDensity,'Nbins',[30 30],'CDataMode','auto','FaceColor','interp');
xlim([0 200E-9]);
ylim([0 100E-9]);
xlabel ('X (m)');
ylabel('Y (m)');
zlabel('Electron Density # of Electrons');
title ('Electron Density Map');















