#create topol

#create topol.top
FILE="prd.gro"
if test  -f $FILE ; then
    echo "$FILE exists"
else 
    echo "$FILE doesn't exist, need to calculate with GROMACS"
  
    cp topol_initial.top topol.top
    echo "create water"
    gmx solvate -cs spc216 -o conf.gro -box 3.152 3.152 3.152 -p topol.top


    

    echo "minimization 1"
    gmx grompp -f mdp/min.mdp -o min -pp min -po min
    gmx mdrun -deffnm min -v

    echo "minimization 2"
    gmx grompp -f mdp/min2.mdp -o min2 -pp min2 -po min2 -c min -t min
    gmx mdrun -deffnm min2 -v

    echo "equilibration 1: temperature coupling"

    gmx grompp -f mdp/eql.mdp -o eql -pp eql -po eql -c min2 -t min2
    gmx mdrun -deffnm eql -v

    echo "equilibration 2: pressure coupling"
    gmx grompp -f mdp/eql2.mdp -o eql2 -pp eql2 -po eql2 -c eql -t eql
    gmx mdrun -deffnm eql2 -v



    echo " production run: NVE ensamble"
    gmx grompp -f mdp/prd.mdp -o prd -pp prd -po prd -c eql2 -t eql2 -p topol.top
    gmx mdrun -deffnm prd -v

fi

#determine number of threads nt to run csg_stat in parallel
nt="$(grep -c processor /proc/cpuinfo 2>/dev/null)" || nt=0
((nt++))

cd Coarse_Graining

echo "Calculating distributions"
csg_stat --top ../prd.tpr --trj ../prd.trr --cg water.xml --options settings.xml --nt $nt --begin 10

echo "Mapping confout.gro to get the starting configuration for coarse-grained runs of ibi/imc/re"
csg_map --top ../prd.tpr --trj ../prd.gro --cg water.xml --out conf.gro 


echo "mapping all the atomistic trajectory to the coarse grained one to perform the optimization of GLE parameters "

csg_map --top ../prd.tpr --trj ../prd.trr --out traj_CG.gro --cg water.xml --vel


echo "mapping the last frame of the atomistic simulation to a CG configuration that will be the initial condition for the GLE simulation"

csg_map --top ../prd.tpr --trj ../prd.trr --out IC.gro --cg water.xml --vel --first-frame 19999 --nframes 1

