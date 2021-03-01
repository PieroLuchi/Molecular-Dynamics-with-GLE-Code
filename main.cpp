#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream> // per scrivere e leggere .txt
#include <cstdlib>
#include <ctime>
#include<stdio.h>

using namespace std;

// VARIABILI GLOBALI

float dt=0.002;
int niter=1750; //1750;
int VAF=1;

float limit=0.8;
int nstout=50;
int nstlist=20;

// termostati

int And=0; //ogni quanti passi applicare andersen

int NH=0; //numero di passi iniziali in cui applicare il Nose_Hoover per termalizzare
float zeta_t=0;
float zeta_hlf_t=0;
float Q=2000;


float bound_up=3.11587;//3.11320;
float bound_down=0;

int num_atom=1001;

float m = 18.0153*0.001; //% [kg/mol] massa CG bead di acqua

int T=298; // K
float kB=0.0083144621*1000;//*0.000001; //  J/(mol K)


float vbolt= 2/sqrt(3.14)*sqrt(2*kB*0.000001*T/m);
float  vxbolt=vbolt/sqrt(3);


// VETTORI
float Ff[500];

float rx[1001];
float ry[1001];
float rz[1001];

float vx[1001];
float vy[1001];
float vz[1001];

float vhlfx[1001];
float vhlfy[1001];
float vhlfz[1001];

int index_vicini[1001];
int num_vicini[1001];
int vicini[140*1001];

float force_x[1001];
float force_y[1001];
float force_z[1001];

/*--LE----*/
int LD=0;
float gammas=0.8;

float par_LE= sqrt(2*kB*0.000001*T*m*gammas*0.86/dt);
float R1;
float R2;
float R3;


/*----GLE-------*/
int GLE=1;
float par_GLE=sqrt(kB*T*0.000001)*1.45;//0.94;//per multi *0.937; //per hybr *0.935;
float par_GLE_k=1;
int M=500;

        float fDx=0;
        float fDy=0;
        float fDz=0;

        float noise_x=0;
        float noise_y=0;
        float noise_z=0;

        float xi_x[500*1001];
        float xi_y[500*1001];
        float xi_z[500*1001];

        float Kx[500];
        float Ky[500];
        float Kz[500];

        float Lx[500];
        float Ly[500];
        float Lz[500];

        float vx_mem[500*1001];
        float vy_mem[500*1001];
        float vz_mem[500*1001];


/*---------------------------------*/



float uniformRandom(){
   //return rand() %10 +1;
  return ( (float)(rand()) + 1. )/( (float)(RAND_MAX) + 1. );
}

float normalRandom(){
    // return a normally distributed random number
  float u1=uniformRandom();
  float u2=uniformRandom();
  return cos(8.*atan(1.)*u2)*sqrt(-2.*log(u1));
}

void load_data(){
ifstream IC("IC.txt");
ifstream force("Ff_ibi_NVE.txt");


    for(int ct=0;ct<num_atom;ct++){
            IC>>rx[ct]>>ry[ct]>>rz[ct]>>vx[ct]>>vy[ct]>>vz[ct];
            vhlfx[ct]=vx[ct];
            vhlfy[ct]=vy[ct];
            vhlfz[ct]=vz[ct];

    }

        for(int ct=0;ct<500;ct++){
            force>>Ff[ct];




    }

    IC.close();
    force.close();

}

void load_data_GLE(){
        ifstream K_file("K.txt");
        ifstream L_file("L.txt");

        for(int ct=0;ct<M;ct++){
            K_file>>Kx[ct]>>Ky[ct]>>Kz[ct];
            L_file>>Lx[ct]>>Ly[ct]>>Lz[ct];
        }

        for(int ct=0;ct<num_atom;ct++){
            xi_x[ct]=normalRandom();//*0.00013971;
            xi_y[ct]=normalRandom();
            xi_z[ct]=normalRandom();

        }


        K_file.close();
        L_file.close();

}

float andersen(){

for(int part=0;part<num_atom;part++){
    vhlfx[part]=normalRandom()*sqrt(kB*T/m*0.86)*0.001; //0.001 è per avere le velocità in nm/ps
    vhlfy[part]=normalRandom()*sqrt(kB*T/m*0.86)*0.001;
    vhlfz[part]=normalRandom()*sqrt(kB*T/m*0.86)*0.001;

    vx[part]=vhlfx[part];
    vy[part]=vhlfy[part];
    vz[part]=vhlfz[part];

}
return 0;
}

void find_neighbors_pbc(){
  int vic=0;
  int neig_count;
  float dx;
  float dy;
  float dz;

float ex=bound_up-limit;

for(int na=0;na<num_atom;na++){
    neig_count=0;
    index_vicini[na]=vic;

    for(int h=0;h<num_atom;h++){
      if(h!=na){
        dx=rx[h]-rx[na];
        dy=ry[h]-ry[na];
        dz=rz[h]-rz[na];

        if ( (abs(dx)<limit || abs(dx)>ex) && (abs(dy)<limit || abs(dy)>ex) && (abs(dz)<limit ||  abs(dz)>ex)){
        //if (sqrt(dx*dx+dy*dy+dz*dz)<limit){

        vicini[vic]=h;
        vic++;
        neig_count++;


        }
        }

    }
    num_vicini[na]=neig_count;

}

}

void find_forces_pbc(){


int at_va;

float dx;
float dy;
float dz;
float distanza;
int dist;
float ex=bound_up-limit;


for(int ct=0;ct<num_atom;ct++){

    force_x[ct]=0;
    force_y[ct]=0;
    force_z[ct]=0;
}

for(int na=0;na<num_atom;na++){



    for(int h=0;h<num_vicini[na];h++){
         at_va=vicini[index_vicini[na]+h];

        if(at_va>na){

            dx=rx[at_va]-rx[na];
            dy=ry[at_va]-ry[na];
            dz=rz[at_va]-rz[na];

            distanza=sqrt(dx*dx+dy*dy+dz*dz);


            if (distanza>ex){

                if (dx>limit+0.1){
                    dx=dx-bound_up;
                }else if (dx<(-limit-0.1)){
                    dx=dx+bound_up;
                }
                if (dy>limit+0.1){
                    dy=dy-bound_up;
                }else if (dy<(-limit-0.1)){
                    dy=dy+bound_up;
                }
                if (dz>limit+0.1){
                    dz=dz-bound_up;
                }else if (dz<(-limit-0.1)){
                    dz=dz+bound_up;
                }

             distanza=sqrt(dx*dx+dy*dy+dz*dz);
            }



            dist=ceil(distanza/0.0020);

            dist=dist-1; // perchè qui la posizione dei vettori inizia da 0!!

            if (dist>500){
                dist=499; // perchè qui la posizione dei vettori inizia da 0!!
            }



                    force_x[na]=force_x[na]+Ff[dist]*dx/distanza;
                    force_y[na]=force_y[na]+Ff[dist]*dy/distanza;
                    force_z[na]=force_z[na]+Ff[dist]*dz/distanza;

                    force_x[at_va]=force_x[at_va]-Ff[dist]*dx/distanza;
                    force_y[at_va]=force_y[at_va]-Ff[dist]*dy/distanza;
                    force_z[at_va]=force_z[at_va]-Ff[dist]*dz/distanza;


                     }

                    }

            }

}

void PBC(int na){
                    if(rx[na]>bound_up){
                        rx[na]=bound_down;
                    }else if( rx[na]<bound_down){
                        rx[na]=bound_up;
                    }


                    if(ry[na]>bound_up){
                        ry[na]=bound_down;
                    }else if( ry[na]<bound_down){
                        ry[na]=bound_up;
                    }

                    if(rz[na]>bound_up){
                        rz[na]=bound_down;
                    }else if( rz[na]<bound_down){
                        rz[na]=bound_up;
                    }
}

void HWBC(int na){

                    if(rx[na]>bound_up){
                        rx[na]=bound_up;
                        vhlfx[na]=-1*vhlfx[na];
                    }else if( rx[na]<bound_down){
                        rx[na]=bound_down;
                        vhlfx[na]=-1*vhlfx[na];
                    }

                    if(ry[na]>bound_up){
                        ry[na]=bound_up;
                        vhlfy[na]=-1*vhlfy[na];
                    }else if( ry[na]<bound_down){
                        ry[na]=bound_down;
                        vhlfy[na]=-1*vhlfy[na];
                    }

                    if(rz[na]>bound_up){
                        rz[na]=bound_up;
                        vhlfz[na]=-1*vhlfz[na];
                    }else if( rz[na]<bound_down){
                        rz[na]=bound_down;
                        vhlfz[na]=-1*vhlfz[na];
                    }
}

float temperatura(float vx[729],float vy[729],float vz[729],float kB, int num_atom,float m){

//float magic = num_atom*vbolt*vbolt*m/2/(kB*0.000001*T)*1000*1000;
float DOF=num_atom*3-3
float somma_vel=0;


    for(int ct=0;ct<num_atom;ct++){
      somma_vel=somma_vel+(vx[ct]*vx[ct]+vy[ct]*vy[ct]+vz[ct]*vz[ct])*1000*1000;
    }



    float Kin=0.5*m*somma_vel;



   float Ti=Kin/(DOF*kB);

   return Ti;

}

void int_zeta(){
float magic =num_atom*vbolt*vbolt*m/2/(kB*T);

float somma_vel=0;
    for(int ct=0;ct<num_atom;ct++){
      somma_vel=somma_vel+(vx[ct]*vx[ct]+vy[ct]*vy[ct]+vz[ct]*vz[ct])*1000*1000;
    }
float Kin=0.5*m*somma_vel;

zeta_hlf_t=zeta_t+dt/(2*Q)*(Kin-magic*kB*T);

somma_vel=0;
    for(int ct=0;ct<num_atom;ct++){
      somma_vel=somma_vel+(vhlfx[ct]*vhlfx[ct]+vhlfy[ct]*vhlfy[ct]+vhlfz[ct]*vhlfz[ct])*1000*1000;
    }
    Kin=0.5*m*somma_vel;
    zeta_t=zeta_hlf_t+dt/(2*Q)*(Kin-magic*kB*T);
}

float sum(float vett[],int leng){
  float somma=0;
    for(int i=0;i<leng;i++){
        somma=somma+vett[i];
    }
        return somma;
}


int main()
{
 /** CARICAMENTO POSIZIONI E VELOCITà INIZIALI E FORZA**/
  load_data();




/** INIZIALIZZO VICINI**/
find_neighbors_pbc();

/** Inizializzo seed generatore casuale*/
srand(time(NULL));
//andersen();
/**Inizializzo velocità**/
//andersen();
cout << temperatura(vhlfx,vhlfy,vhlfz,kB,num_atom,m)<<endl;


/** GLE**/
        if(GLE==1){
          load_data_GLE();
          And=0;
         cout << par_GLE <<endl;
cout << Kx[0] << " "<<Ky[0]<<" " <<Kz[0]<<endl;
cout << Lx[0] << " "<<Ly[0]<<" " <<Lz[0]<<endl;
        }


//        find_forces_pbc();
//
//        for(int ct=0;ct<num_atom;ct++){
//         cout << force_x[ct]<<force_y[ct]<<force_z[ct]<<endl;
//        }

/**INIZIALIZZAZIONE FILE SALVATAGGIO**/

    remove("gro.txt");
    ofstream gro("gro.txt");
    gro << "niter" << " "<< "nstout" << " " <<"dt" << " " << "num_atom" << " " << "bound_up"<<endl;
    gro << niter << " "<< nstout << " " <<dt << " " << num_atom << " " << bound_up<<endl;
    gro << "x     y     z"<<endl;

    remove("veloc.txt");
    ofstream veloc("veloc.txt");
    veloc << "vx     vy     vz"<<endl;


    for(int z=0;z<num_atom;z++){
        gro << rx[z]<< " " << ry[z]<<" "<<rz[z]<<endl;
        veloc << vx[z]<< " " << vx[z]<<" "<<vx[z]<<endl;
    }

    remove("Temp_data.txt");
    ofstream Temp_data("Temp_data.txt");

    Temp_data<< " Temperatura"<<endl;


/**INTEGRAZIONE EQUAZIONI DEL MOTO*/

    /*contatori vari*/


    float inv_dt=1.0/dt;
    float inv_m=1.0/m;

        int count_nstlist=1;
        int count_temp=1;
        int count_nstout=1;

        float Ti;
        int volte_qui=0;
        float  Temp_media=0;


        cout<<"begin"<<endl;


    for(int t=1;t<=niter;t++){

        if(VAF==1){
            if(t<1001){
                nstout=1000;
            }

            if(t==1001){
                nstout=1;
            }
        }

            /**Termostato di Andersen**/

            if(count_temp==And){
                andersen();
                count_temp=0;
                cout<<" reset temp" <<endl;
            }


            /** Aggiorna Vicini*/
                  if (count_nstlist== nstlist){
                      find_neighbors_pbc();
                      count_nstlist=0;
                  }


             /**Calcolo Forze**/

                 find_forces_pbc();



                 /** Termostato Nose-Hoover**/
                 if(t<NH ){
                 int_zeta();
                 }

                 if(t==NH){
                 zeta_t=0;
                 }


        /**Velocity Verlet Integrator**/
if (GLE==0){

        if(LD==0){

                float den=1/(1+dt*0.5*zeta_t);

           for(int na=0;na<num_atom;na++){


                    if(t>1){

                    vx[na]=(vhlfx[na]+0.5*force_x[na]*dt*inv_m)*den;
                    vy[na]=(vhlfy[na]+0.5*force_y[na]*dt*inv_m)*den;
                    vz[na]=(vhlfz[na]+0.5*force_z[na]*dt*inv_m)*den;
                    }

                 vhlfx[na]=vx[na]+0.5*force_x[na]*dt*inv_m-zeta_t*vx[na]*dt*0.5;
                 vhlfy[na]=vy[na]+0.5*force_y[na]*dt*inv_m-zeta_t*vy[na]*dt*0.5;
                 vhlfz[na]=vz[na]+0.5*force_z[na]*dt*inv_m-zeta_t*vz[na]*dt*0.5;

                 rx[na]=rx[na]+vhlfx[na]*dt;
                 ry[na]=ry[na]+vhlfy[na]*dt;
                 rz[na]=rz[na]+vhlfz[na]*dt;

                       PBC(na);
                    }
        }

       if(LD==1){
                And=0;
                float den=1/(1+dt*0.5*gammas);

         for(int na=0;na<num_atom;na++){


                 R1=normalRandom();
                 R2=normalRandom();
                 R3=normalRandom();

                        if(t>1){

                    vx[na]=(vhlfx[na]+0.5*(force_x[na]+par_LE*R1)*dt*inv_m)*den;
                    vy[na]=(vhlfy[na]+0.5*(force_y[na]+par_LE*R2)*dt*inv_m)*den;
                    vz[na]=(vhlfz[na]+0.5*(force_z[na]+par_LE*R3)*dt*inv_m)*den;
                    }

                   vhlfx[na]=vx[na]+0.5*(force_x[na]-gammas*m*vx[na]+par_LE*R1)*dt*inv_m;
                   vhlfy[na]=vy[na]+0.5*(force_y[na]-gammas*m*vy[na]+par_LE*R2)*dt*inv_m;
                   vhlfz[na]=vz[na]+0.5*(force_z[na]-gammas*m*vz[na]+par_LE*R3)*dt*inv_m;


                    rx[na]=rx[na]+vhlfx[na]*dt;
                    ry[na]=ry[na]+vhlfy[na]*dt;
                    rz[na]=rz[na]+vhlfz[na]*dt;

                       PBC(na);
                 }
              }

}else if (GLE==1){

    float denx=1/(1+Kx[0]*dt*dt*0.25*inv_m);
    float deny=1/(1+Ky[0]*dt*dt*0.25*inv_m);
    float denz=1/(1+Kz[0]*dt*dt*0.25*inv_m);

            for(int na=0;na<num_atom;na++){


                        /**Memory Kernel e Colored Noise**/

                          fDx=0;
                          fDy=0;
                          fDz=0;
                          noise_x=0;
                          noise_y=0;
                          noise_z=0;

                                for(int jj=1; jj<min(M,t-1);jj++){
                                     int ll=(jj-1)*num_atom;

                                        fDx=fDx+dt*Kx[jj]*vx_mem[na+ll];
                                        fDy=fDy+dt*Ky[jj]*vy_mem[na+ll];
                                        fDz=fDz+dt*Kz[jj]*vz_mem[na+ll];

                                        noise_x=noise_x+Lx[jj]*xi_x[na+ll];
                                        noise_y=noise_y+Ly[jj]*xi_y[na+ll];
                                        noise_z=noise_z+Lz[jj]*xi_z[na+ll];

                                }


                    if(t>1){

                    vx[na]=(vhlfx[na]+0.5*dt*inv_m*(force_x[na]-fDx+par_GLE*noise_x))*denx;
                    vy[na]=(vhlfy[na]+0.5*dt*inv_m*(force_y[na]-fDy+par_GLE*noise_y))*deny;
                    vz[na]=(vhlfz[na]+0.5*dt*inv_m*(force_z[na]-fDz+par_GLE*noise_z))*denz;
                    }

                 vhlfx[na]=vx[na]+0.5*dt*inv_m*(force_x[na]-0.5*Kx[0]*vx[na]*dt-fDx+par_GLE*noise_x);
                 vhlfy[na]=vy[na]+0.5*dt*inv_m*(force_y[na]-0.5*Ky[0]*vy[na]*dt-fDy+par_GLE*noise_y);
                 vhlfz[na]=vz[na]+0.5*dt*inv_m*(force_z[na]-0.5*Kz[0]*vz[na]*dt-fDz+par_GLE*noise_z);

                 rx[na]=rx[na]+vhlfx[na]*dt;
                 ry[na]=ry[na]+vhlfy[na]*dt;
                 rz[na]=rz[na]+vhlfz[na]*dt;

                       PBC(na);

           }

}

               /**Storage velocity for GLE**/
              if(GLE==1){


                    for(int ct=(num_atom*M)-1; ct>=num_atom;ct--){

                       vx_mem[ct]=vx_mem[ct-num_atom];
                       vy_mem[ct]=vy_mem[ct-num_atom];
                       vz_mem[ct]=vz_mem[ct-num_atom];

                       xi_x[ct]=xi_x[ct-num_atom];
                       xi_y[ct]=xi_y[ct-num_atom];
                       xi_z[ct]=xi_z[ct-num_atom];

                    }

                    for(int ct=0;ct<num_atom;ct++){
                       vx_mem[ct]=vhlfx[ct];
                       vy_mem[ct]=vhlfy[ct];
                       vz_mem[ct]=vhlfz[ct];

                       xi_x[ct]=normalRandom();
                       xi_y[ct]=normalRandom();
                       xi_z[ct]=normalRandom();

                    }

//                    for(int g=0;g<num_atom*2;g++){
//
//                       // cout<<vx_mem[g]<<endl;
//                       vel_mem<<vx_mem[g]<<endl;
//
//                       if (g==num_atom-1){
//                      vel_mem<<"------------------------------------"<<endl;
//                    }
//
//                    }
//                      vel_mem<<"------------------------------------"<<endl;
//                    vel_mem<<"------------------------------------"<<endl;
//
//
             }



    /** Write data on txt file and display informations**/
      if (count_nstout==nstout){

            cout<<"Time step="<<t<<endl;

          for(int z=0;z<num_atom;z++){
                gro << rx[z]<< " " << ry[z]<<" "<<rz[z]<<endl;
                veloc << vx[z]<< " " << vy[z]<<" "<<vz[z]<<endl;
            }



            Ti=temperatura(vhlfx,vhlfy,vhlfz,kB,num_atom,m);
            cout<<"Temperatura istantanea=" << Ti<<endl;
            Temp_data << Ti<<endl;


//          volte_qui=volte_qui+1;
          Temp_media=Temp_media+Ti;
//
//           cout<< " Temperatura media ="<<Temp_media/volte_qui<<endl;
//


          count_nstout=0;

       }


       /*Aggiornamento contatori*/
         count_temp=count_temp+1;
         count_nstlist=count_nstlist+1;
         count_nstout=count_nstout+1;


    }

cout<< "TEMPERATURA MEDIA="<<Temp_media/niter*nstout;
cout << "FINE"<<endl;
}
