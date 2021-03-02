clearvars -except L1
close all

%funcform-------
% load('C:\Users\piero\Dropbox\AA Collaborazione\Memory Kernel Opt\Real Kernel\Kottim')
load('C:\Users\piero\Dropbox\AA Collaborazione\Memory Kernel Opt\Memory Kernel water functional form\K_ottim_func_form')
%free par-------
% load('C:\Users\piero\Dropbox\AA Collaborazione\Memory Kernel Opt\Memory Kernel water free param method\K_ottim_free_par')
% load('C:\Users\piero\Dropbox\AA Collaborazione\Memory Kernel Opt\Memory Kernel water free param method\L_ottim_free_par')
existL=0;
%---------------

fake=1;
backward=0;

dt=0.008
salta=dt*1000/2
M=length(K)
T=298; %K
kB=0.0083144621*1000*10^-6;  %J/(mol K)

stop_cond=0.006; 
diff_cond=0.05;

if M==150
stop_cond=0.06; %0.006 per M=350
diff_cond=0.5;
end

%% processo inverso tolgo valori per ottimizzare più velocemente
if backward==1
    M=M/salta;
    count=1;
    for i=1:salta:M*salta
        Kn(count,1:3)=K(i,:);
        count=count+1;
    end    
    K=Kn;
    
end

%% se ho già i parametri L

if existL==1
  
L1=L(:,1);
L2=L(:,2);
L3=L(:,3);  
    
end

%% ottimizzazione
if existL==0
[L1]=Opt_Noise_par(K(1:M,1),M,1,stop_cond,diff_cond);
fprintf('L1 fatto \n')
[L2]=Opt_Noise_par(K(1:M,2),M,L1,stop_cond,diff_cond);
fprintf('L2 fatto \n')
[L3]=Opt_Noise_par(K(1:M,3),M,L1,stop_cond,diff_cond);
% L2=L1;
% L3=L1;
% [L1]=Opt_Noise_par(K(1:M,1),M,L1,stop_cond);
% [L2]=Opt_Noise_par(K(1:M,2),M,L2,stop_cond);
% [L3]=Opt_Noise_par(K(1:M,3),M,L3,stop_cond);
end
%%
if fake==1
    count=1;
   L=zeros(M*salta,3);
   Kn=zeros(M*salta,3);
   
   if dt==0.004
   magic=0.71 %0.004
   elseif dt==0.006
   magic=0.575 %0.006
   elseif dt==0.008
   magic=0.5 %0.008
   end
   
   for m=1:salta:M*salta
   L(m,1)=L1(count)*magic;
   L(m,2)=L2(count)*magic;
   L(m,3)=L3(count)*magic;
   
   Kn(m,1)=K(count,1);
   Kn(m,2)=K(count,2);
   Kn(m,3)=K(count,3);
   count=count+1;
   end
   
   if salta==2
   for m=2:salta:salta*M-1
       
       L(m,1)=(L(m-1,1)+L(m+1,1))/2;
       L(m,2)=(L(m-1,2)+L(m+1,2))/2;
       L(m,3)=(L(m-1,3)+L(m+1,3))/2;
       
       Kn(m,1)=(Kn(m-1,1)+Kn(m+1,1))/2;
       Kn(m,2)=(Kn(m-1,2)+Kn(m+1,2))/2;
       Kn(m,3)=(Kn(m-1,3)+Kn(m+1,3))/2;
   end
   end
   
   if salta==3
       
       for m=2:salta:salta*M-2
       L(m,1)=(L(m-1,1)+L(m+2,1))/2;
       L(m,2)=(L(m-1,2)+L(m+2,2))/2;
       L(m,3)=(L(m-1,3)+L(m+2,3))/2;
       
       Kn(m,1)=(Kn(m-1,1)+Kn(m+2,1))/2;
       Kn(m,2)=(Kn(m-1,2)+Kn(m+2,2))/2;
       Kn(m,3)=(Kn(m-1,3)+Kn(m+2,3))/2;
       end
       
       for m=3:salta:salta*M-1
       L(m,1)=(L(m-1,1)+L(m+1,1))/2;
       L(m,2)=(L(m-1,2)+L(m+1,2))/2;
       L(m,3)=(L(m-1,3)+L(m+1,3))/2;
       
       Kn(m,1)=(Kn(m-1,1)+Kn(m+1,1))/2;
       Kn(m,2)=(Kn(m-1,2)+Kn(m+1,2))/2;
       Kn(m,3)=(Kn(m-1,3)+Kn(m+1,3))/2;
       end
   end
   
   if salta==4
       for m=3:salta:salta*M-2
           
           for p=1:3
       L(m,p)=(L(m-2,p)+L(m+2,p))/2;
       Kn(m,p)=(Kn(m-2,p)+Kn(m+2,p))/2;
           end

       end
       
       for p=1:3
       L(M*salta-1,p)=(L(M*salta-3,p))/2;
       Kn(M*salta-1,p)=(Kn(M*salta-3,p))/2;
       end
       
       for m=2:salta:salta*M-1
          for p=1:3
       L(m,p)=(L(m-1,p)+L(m+1,p))/2;
       Kn(m,p)=(Kn(m-1,p)+Kn(m+1,p))/2;
       end
       
       end
       
       for m=4:salta:salta*M-1
           for p=1:p
       L(m,p)=(L(m-1,p)+L(m+1,p))/2;      
       Kn(m,p)=(Kn(m-1,p)+Kn(m+1,p))/2;
           end
       end
       
       for p=1:p
       L(M*salta,p)=(L(M*salta-1,p))/2;      
       Kn(M*salta,p)=(Kn(M*salta-1,p))/2;
           end
       
       
   end
   
   
end

if fake==0
    L=[L1 L2 L3]
    Kn=K
end
%%
f1=xcorr(L(:,3));
% plot(0:0.004:0.004*(175-1),f1(M:end))
plot(0:0.002:0.002*(salta*M-1),f1(salta*M:end))
hold on
plot(0:0.002:0.002*(salta*M-1),Kn(:,3))
legend('corr','K')

%% 
%  save('L_fin','L')
  K=Kn;
  save('K_ottim_extended','K')
  save('L_ottim_extended','L')
 %%
fprintf('salvati in MD \n')
fileID = fopen('D:\Git Repositories\Molecular-Dynamics-with-GLE-Code\K.txt','w');
fprintf(fileID,'%.16f %24.16f %24.16f\n',K');
fclose(fileID);

fileID = fopen('D:\Git Repositories\Molecular-Dynamics-with-GLE-Code\L.txt','w');
fprintf(fileID,'%.16f %24.16f %24.16f\n',L');
fclose(fileID);




 