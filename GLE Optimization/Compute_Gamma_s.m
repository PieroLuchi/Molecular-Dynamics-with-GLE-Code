function Gamma_s=Compute_Gamma_s(num_atom,M,N,traj)
Gamma_s=zeros(3,M);
 for atomo=1:6:num_atom*3
     
    (atomo+5)/6
        vx=traj(1:N,atomo+2+1);
        vy=traj(1:N,atomo+2+2);
        vz=traj(1:N,atomo+2+3);
                for n=M:N
                    Gamma_s(1,:)=Gamma_s(1,:)+pseudoxcorr(vx,vx,n,M);
                    Gamma_s(2,:)=Gamma_s(2,:)+pseudoxcorr(vy,vy,n,M);
                    Gamma_s(3,:)=Gamma_s(3,:)+pseudoxcorr(vz,vz,n,M);
                end
     
   
 end 