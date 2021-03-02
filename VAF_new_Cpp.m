%% Velocity Autocorrelation

% ricordarsi di usare una traiettoria molto fine !!! 750 step a 0.002 ps
    close all
    clear all
 [niter nstout dt num_atom bound_up]=textread("gro.txt",' %f %f %f %f %f ',1,'headerlines',1);

L=bound_up;

%% carico dati
%scelgo quanti frame iniziali salatare  
  f=num_atom*250;
  g=num_atom*1000;
  
   head=1;
   [x y z]=textread("veloc.txt",' %f %f %f ',g,'headerlines',head);
   v=[x,y,z];
   


%calcolo quindi il numero di frame
  frames=length(v)/num_atom

 %modifico topologia dati
  vv=zeros(num_atom,3,frames);
  
  for t=1:frames
      vv(:,:,t)=v(1+(t-1)*num_atom:t*num_atom,:);
  
  end
  
%% calcolo VAF 3D

tempo_max=1.5/(dt*1); %calcolo il VAF solo fino a 1.5 ps
VAF=zeros(tempo_max,3);
Norm=[0 0 0];
for t0=1:frames
  for t=1:min(frames-t0+1,tempo_max)
      
      VAF(t,:)=VAF(t,:)+sum(vv(:,:,(t-1)+t0).*vv(:,:,t0));
    
  end
  Norm=Norm+sum(vv(:,:,t0).*vv(:,:,t0));
end
 


  VAF=VAF./Norm;%sum(vv(:,:,t0).*vv(:,:,t0));
  
  VAF=sum(VAF')/3;
  
  DIFF_COEFF=sum(VAF*dt)*10^2
  
  
   %% Dati Termodinamici 

T=textread("Temp_data.txt",' %f ',1000,'headerlines',2);




TEMP_MEDIA=mean(T(10:end))

ERRORE_TEMP = std(T)/sqrt(length(T)) 
  
  
  %% grafici
%  load('C:\Users\piero\Dropbox\AA Collaborazione\Memory Kernel Opt\Memory Kernel water functional form\K_ottim_func_form')
% linspace(0,750,100),K(:,1)/(K(1,1)) 
load('VAF_atomistic')
  figure
  plot(0:dt*nstout:(tempo_max-1)*dt*nstout,[VAF' VAF_atomistic])
  hold on
  %plot(0:dt*nstout:(400-1)*dt*nstout,[Kx./Kx(1)])
  title('3D velocity Autocorr Function')
  xlabel('time [ps]')
  %% salva VAF
%save('C:\Users\piero\Dropbox\AA Collaborazione\Codice Ottimizz Kernel TNL Finale\Prod Grafici\VAF_s2','VAF')
 
  %[Kx Ky Kz]=textread("K.txt",' %f %f %f ',400,'headerlines',0);

  