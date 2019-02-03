function[electrons] = move3Final(electrons)
global dt Pscat NumP Vth  diffusive boxes specular ax ax1
    boxes;
    a=zeros(NumP,1);
    
    top= electrons(:,2)+electrons(:,4)*dt*2 >100E-9; % logical indexing for if the electron goes above the top limit
    bottom =electrons(:,2)+electrons(:,4)*dt*2<0; %logical indexing for if the electron goes below the bottom limit
    
    sp=rand(NumP,1);
    Scatter = Pscat(:,1) >sp;
    
    [a(:,1),ax]=logixbox(electrons);
  
 
  random1 = rand()*10;
  if random1 >5
      random1= random1/5;
  end
  
  random2 = rand()*10;
  if random2>5
      random2=random2/5;
  end
   
   
    electrons(:,4)=electrons(:,4).*~Scatter ... %  Vy = Vy unless scattering
            - 2*electrons(:,4).*top.*specular -2*electrons(:,4).*bottom.*specular ... % Specular reflection condition top and bottom
             - random1*electrons(:,4).*top.*diffusive -random2*electrons(:,4).*bottom.*diffusive ... %Diffusive reflection condition top and bottom
            + randn()*(Vth/sqrt(2)).*~top.*~bottom.*Scatter.*~a(:,1) ... % Scattering condition
            -2*electrons(:,4).*a(:,1).*~ax(:,1).*specular ... % Specular box conditions
            -random1*electrons(:,4).*a(:,1).*~ax(:,1).*diffusive ; % Diffusive box conditions
    
    
    % Logical indexing for top and bottoms and also scattering condition
    
    electrons(:,3)=electrons(:,3).*~Scatter ... % Vx=Vx unless scattering
            + randn()*(Vth/sqrt(2)).*~top.*~bottom.*Scatter.*~a(:,1) ... % scattering condition
            -2*electrons(:,3).*a(:,1).*ax(:,1).*specular... % Specular box conditions
            - random2*electrons(:,3).*a(:,1).*ax(:,1).*diffusive ; %diffusive box conditions
  
    %-2*electrons(:,3).*a*specular.*~Scatter -rand()*electrons(:,3).a*diffusion.*~Scatter
    % Scattering condution for x 
    
    electrons(:,1)=mod(electrons(:,1)+electrons(:,3)*dt,200E-9); 
    % X position update and wrap around condition for cts boundaries
    
    
    electrons(:,2)=electrons(:,2)+electrons(:,4)*dt;
    % Y position update
    
end


