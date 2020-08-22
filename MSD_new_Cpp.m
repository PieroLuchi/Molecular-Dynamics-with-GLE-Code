clear all
close all

[niter nstout dt num_atom bound_up]=textread("gro.txt",' %f %f %f %f %f ',1,'headerlines',1);

L=bound_up;


%% carico dati posizioni

%scelgo quanti frame iniziali salatare   
  f=num_atom*0; %inizio
  g=num_atom*750;  %quanti frame caricare dal file e analizzare

   head=3+f;
%    [x y z]=textread("gro.txt",' %f %f %f ',(niter/nstout)*num_atom,'headerlines',head);
[x y z]=textread("gro.txt",' %f %f %f ',g,'headerlines',head);
   r=[x,y,z];
   

% calcolo quindi il numero di frame orrisponenti
  frames=length(r)/num_atom
  
  % metto i frame in un vettore 3d più comodo per i calcoli
  rr=zeros(num_atom,3,frames);
  
  for t=1:frames
      rr(:,:,t)=r(1+(t-1)*num_atom:t*num_atom,:);
  
  end
  
 %definisco vettori per grafico 
  MSDt=0:dt*nstout:(frames-1)*nstout*dt'; %ps
  MSDt=MSDt';
  MSD=zeros(length(MSDt),1);
  
  
  %% per le PBC rettifico le traiettorie:
  for double=1:10
for t=1:length(MSD)-1
    
       differ=rr(:,:,t+1)-rr(:,:,t);
    
          for i=1:num_atom
             for comp=1:3
                if differ(i,comp)>L/2
                    rr(i,comp,t+1)=rr(i,comp,t+1)-L;
                elseif differ(i,comp)<(-L/2)
                    rr(i,comp,t+1)=rr(i,comp,t+1)+L;
                end
             end
            end

end
  end
% controllo traiettoria
figure
for t=1:length(MSD)
plot(t,rr(200,1,t),'.');
hold on
end

%% calcolo diffusion coefficient
for t=1:length(MSD)
    
   val(t,1:3)= (sum((rr(:,:,t)-rr(:,:,1)).^2)/num_atom);
    
   
end
MSD=sum(val')';

fit_num=5;floor(frames/2); % inizio fitting
[MSDpar,b]=polyfit(MSDt(fit_num:end),MSD(fit_num:end),1);

figure
plot(MSDt(1:end),MSD(1:end),'.',MSDt(fit_num:end),polyval(MSDpar,MSDt(fit_num:end)))

DIFF_COEFF=MSDpar(1)/6*10^(-2)*10^5 %cm^2/s

%std(diff(MSD)./diff(MSDt))/sqrt(length(diff(MSD)))

ERR_DIF_COEFF=std((MSD(fit_num:end)-polyval(MSDpar,MSDt(fit_num:end)))/6*10^(3))/sqrt(length(MSD))

%% Dati Termodinamici 

T=textread("Temp_data.txt",' %f ',niter/nstout,'headerlines',2);




TEMP_MEDIA=mean(T(10:end))

ERRORE_TEMP = std(T)/sqrt(length(T))
  
