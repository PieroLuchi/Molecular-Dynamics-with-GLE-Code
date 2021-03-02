%% LOAD DATA from GROMACS output file of the trajectories  
% this function takes the .gro format Gromacs input dile of CG trajectories
% and save them in a simpler format that is used in our GLE optimization

clc 
clear all
close all

% parameters
L=3.1132 %lenght of box 
frames=15000 %% frames of the trajectory
num_atom=1001*3  % number of atom in atomistic config
num_atom_CG=1001 % number of atom in CG config

%position of the Gromacs file of trajectories in .gro format 
Gromacs_traj='C:\Users\piero\Desktop\traj\traj_60000_008_NVE.txt'

%% load file
count=0;

tic
head=0
[x y z vx vy vz]=textread(Gromacs_traj,'%*s %*s %*s %f %f %f %f %f %f',num_atom_CG*frames,'headerlines',head);
coords=[x,y,z,vx,vy,vz];
toc

% save('coords','coords')
%%
% load('coords')

%% Clean and prepare trajectories

tic
atom=1;
traj=zeros(frames,num_atom_CG*6);
for i=1:6:num_atom_CG*6
    
    fr=0;
        for j=atom:num_atom_CG:num_atom_CG*frames
            fr=fr+1;
            traj(fr,i:i+5)=coords(j,:);
        end
        atom=atom+1;
end
toc

 plot(traj(:,105))
hold on


%save('traj','traj')


%% Rettify trajectories
% Since they are computed in a periodic boundary condition framework,
% We need to "rectify" these trajectory to avoid non-pysical jump form 
% a face to the other of the simulation box
 
%load('traj')

%%
tic
for red=1:3
for na=1:6:num_atom_CG*6
 for t=1:frames-1

    diff=traj(t+1,na:na+2)-traj(t,na:na+2);
    for comp=1:3
                if diff(comp)>L-0.5
                    traj(t+1,na-1+comp)=traj(t+1,na-1+comp)-L;
                elseif diff(comp)<(-L+0.5)
                    traj(t+1,na-1+comp)=traj(t+1,na-1+comp)+L;
                end
    end
 end
end
end
toc

 plot(traj(:,105))
 
 save('traj','traj')
 
%% Separation of velocities and positions
% uncomment this section if you want to create separate files 
% for velocity and position.  Not needed in the present analysis


% traj_vel=zeros(frames,num_atom_CG*3);
% traj_pos=zeros(frames,num_atom_CG*3);
% count=1
% for na=1:6:num_atom_CG*6
%     
%     traj_pos(:,count:count+2)=traj(:,na:na+2);
%     traj_vel(:,count:count+2)=traj(:,na+3:na+5);
%     count=count+3;
% end
% 
% save('traj_pos','traj_pos');
% save('traj_vel','traj_vel');
%  
%  %%
% for i=1:10000
%  plot3(traj_pos(:,1),traj_pos(:,2),traj_pos(:,3))
%  axis([0,L+1,0,L+1,0,L+1])
%  grid on
% % 
% end
