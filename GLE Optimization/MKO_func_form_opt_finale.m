%% Kernel Optimization 
%this is the main code. This code takes the trajectories obtaine 
% in "Load_trajectories.m" and perform the optimization procedure
% described in the work

clearvars -except traj
close all


%takes the trajectories

%addpath('C:\Users\piero\Dropbox\AA Collaborazione\Memory Kernel Opt')
 %load('D:\Documenti\DOTTORATO\Traiettoria per ottimizzazione parametri\traj_60000_008_NVE')
%% variables
%% Attention: this are the values that should be fixed if you are using another sistem!
dt=0.008 %time step in ps
bound_up=3.0932; % side of the simulation box in nm
num_atom=1001;
m = 18.0153*0.001; % kg/mol massa CG bead of water
T=298; %temperature in K
% kB=0.0083144621*1000;%*10^-6 *m %J/(mol K)
% kB=kB*10^(-6);

kB=8.3144621*10^-6 % kb in (J/(mol K)) tranformed 
                    % in kg/mol nm^2/ps^2 for the code

D=2.44*10^(-3)% Diffusion coefficient in nm^2/ps

M=150 % memory kernel depth
N=15000 % number of timestep to take


%% Compute Gamma_s

% compute the Gamma described in work. modifiy the name in necessary.
% if Gamma_s is already computed it does not compute it again
existence_of_Gamma_s_mat=exist('Gamma_s_150.mat')
existence_of_Gamma_s=exist('Gamma_s_150')

if (existence_of_Gamma_s==0 && existence_of_Gamma_s_mat==0)
    %the variable does not exists, compute Gamma_s
    tic  
    Gamma_s=Compute_Gamma_s(num_atom,M,N,traj)
    toc
    save('Gamma_s_150','Gamma_s')
elseif existence_of_Gamma_s==1
    %the variable exists in workspace, do nothing
elseif (existence_of_Gamma_s==0 && existence_of_Gamma_s_mat==2)
    %load already computed Gamma_S 
    load('Gamma_s_150')
    
end

%% OPTIMIZATION
%Two stages optimizatio. 
% First stage: only functional form g(t)
% Second stage: fix functional form with optimizeed local
%               parameters (b(t) function)
%               
D_par=10
%% STAGE 1
figure
% prepare parameters vector
par_s1=zeros(3,2);

%time
t=0:dt:dt*(M-1);

%functional form of the stage one
K=@(P) P(1)*exp(-P(2)*t.^2);

%loop over the three spatial directions, k= x,y,z = spatial axes    
for k=1:3
    % optimization needed to make the kernel to cover all the
    % time dt*M it covers
    for g=1:10
        
    %Define the lagrangian action
    Sigma=@(P) -dt^3*Gamma_s(k,:)*K([P])';

    %Minimization parameters:--------
        % constraint-----------------
        elem=kB*T/((D)*D_par);
        %equality constraint ceq=0
        ceq=@(P) [sum(K([P])*dt)-elem];
        %inequality constraint c<=0
        c=@(P) 0

        nonlincon= @(P)deal(c(P),ceq(P));
        %----------------------------
        %fmincon parameters:
        x0 = [ 0  1 ];%initial guess
        A = []; % No ther constraints
        b = []; 
        Aeq = [];
        beq = [];
        lb =  [ 0 0  ]; %lower bound for the values
        ub = [1.1-0.1*g 100]; %upper bound for the values
        %---------------------
    % minimization:
       par_s1(k,:)=fmincon(Sigma,x0,A,b,Aeq,beq,lb,ub,nonlincon)
       Kcontrol=K(par_s1(k,1:end));
        if Kcontrol(end-40)>0.01
            break
        end
    
    end
end

%Plot of the kernel
plot(t,K(par_s1(1,:)),'.-')
legend('stage1')
xlabel('time [ps]')
hold on 
pause(0.2)
% 

Kottim=[K(par_s1(1,1:end))' K(par_s1(2,1:end))' K(par_s1(3,1:end))'];

L1_norm=(sum(K(par_s1(1,1:end))));
L2_norm=sqrt(sum(K(par_s1(1,1:end).^2)));

% L1_norm_dt=(sum(K(par_s1(1,1:end))*dt))
% L2_norm_dt=sqrt(sum(K(par_s1(1,1:end).^2)*dt))
% %L3_norm=(sum(K(par_s1(1,1:end).^3)))^(1/3)

%save the shape of the first stage Kernel if needed
% K_s1=K(par_s1(1,1:end))
% save('K_s1','K_s1')
%% STAGE 2  
%prepare the vector for stage 2 parameters
par_s2=zeros(3,M+2);

lambda=L1_norm %initial lambda



%loop over the three spatial directions, k= x,y,z = spatial axes    
for k=1:3  
    
    % Memory Kernel functional form of the stage two
    K=@(P) P(1)*exp(-P(2)*t.^2).*P(3:M+2);
    % Lagrangian Action
    Sigma=@(P) -dt^3*Gamma_s(k,:)*K([P])';  
    %lambda=lambda-0.2
    

    %Minimization parameters:---------------
        % constraint------------------------
        elem=kB*T/((D)*D_par)
        ceq=@(P)[sum(K([P])*dt)-elem]

        c=@(P) [sqrt(sum((K([P])).^2))-lambda;...
       sqrt(sum(P(3:M+2).^2))-(L1_norm/L2_norm)];

       %c = @(P) 0;
        nonlincon= @(P)deal(c(P),ceq(P));
        %-------------------------------------
        %fmincon parameters:
        x0 = [ par_s1(k,:) zeros(1,M)-1];
        A = []; % No other constraints
        b = [];
        Aeq = [];
        beq = [];
        %lb = [ par_s1(k,:) -(1+20*exp(-5*linspace(0,2,M)))  ];
        %ub = [ par_s1(k,:) (1+20*exp(-5*linspace(0,2,M))) ]
        %lb = [ par_s1(k,:) -1*(10*exp(-0.4*linspace(1,15,M)))  ];
        %ub = [ par_s1(k,:)  2(10*exp(-0.4*linspace(1,15,M))) ];
        lb = [ par_s1(k,:) zeros(1,M)-1 ];
        ub = [ par_s1(k,:)  zeros(1,M)+1000];
        options = optimoptions('fmincon','MaxFunctionEvaluations',10000,'MaxIterations',100000);

    % minimization:    
        par_s2(k,:)=fmincon(Sigma,x0,A,b,Aeq,beq,lb,ub,nonlincon,options);


end


b=par_s2(1,3:end);
b_L2=sum(b.^2)^(1/2)
L2=sqrt(sum(K(par_s2(1,:)).^2))
%Kottim=[Km',Km',Km'];

plot(t,K(par_s2(1,:)),'.-')
hold on
legend('stage1','stage2')
xlabel('time [ps]')
hold on 
pause(0.2)


Kottim=[K(par_s2(1,:))' K(par_s2(1,:))' K(par_s2(1,:))'];



% save data
% K_s2=K(par_s2(1,1:end))
% save('K_s2','K_s2')
%% SAVE KERNEL
% here the value of the K that are zero are removed to
% help the optimization of noise parameters. 
% then we save it if necessary  

if Kottim(125)<=0.04
    M=125
    K=Kottim(1:M,:);
elseif Kottim(100)<=0.01
    M=100
    K=Kottim(1:M,:);

end
save('K_ottim','K')

%% Find L: Noise parametes

% here we find the noise parameters with the optimization procedure 
% described in the paper

%% Data
load('K_ottim')
vacancies=1;% =1 if the code should add, in an appropriate way, some point to the Kenrnel           % to 
            % to make the kernel suitable for simulation with delta t of
            % 0.002 ps 
salta=dt*1000/2 % point that lack in the kernel ( with dt=0.008 the K lacks 4
                % 3 points between each couple of points
                
stop_cond=0.006; % stop condition: stop minimization if the error is less than this
diff_cond=0.05; % stop condition to avoid oscillatory non convergent resutls
%                : Stop the optimization if the difference between
                 %two subsequent optimization errors 
                 %is less than this



%% Optimization 
% L1,L2,L3 are the colored noise parameters
[L1]=Opt_Noise_par(K(1:M,1),M,1,stop_cond,diff_cond);
fprintf('L1 done \n')
% I use L1 as initial guess here to faster convergence
% [L2]=Opt_Noise_par(K(1:M,2),M,L1,stop_cond,diff_cond);
% fprintf('L2 fatto \n')
% [L3]=Opt_Noise_par(K(1:M,3),M,L1,stop_cond,diff_cond);
% fprintf('L3 fatto \n')
L2=L1
L3=L1

%% Enhance Kenrel and Noise

%here a Kernel and noise for a fine 0.002 ns simulation is constructed
if vacancies==1
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

if vacancies==0
    L=[L1 L2 L3]
    Kn=K
end
%% Plot the Kernel and the optimized noise parameter L
close all
f1=xcorr(L(:,1));
% plot(0:0.004:0.004*(175-1),f1(M:end))
plot(0:0.002:0.002*(salta*M-1),f1(salta*M:end))
hold on
plot(0:0.002:0.002*(salta*M-1),Kn(:,1))
legend('corr','K')


%%  Save data
%save in the folder the .mat file for L and K
K=Kn;
save('K_ottim_extended','K')
save('L_ottim_extended','L')
%%
%save in a coustomized directory the txt file for K and L to be used in the
%simulation with the GLE solver
fprintf('salvati in MD \n')
fileID = fopen('..\K.txt','w');
fprintf(fileID,'%.16f %24.16f %24.16f\n',K');
fclose(fileID);

fileID = fopen('..\L.txt','w');
fprintf(fileID,'%.16f %24.16f %24.16f\n',L');
fclose(fileID);
