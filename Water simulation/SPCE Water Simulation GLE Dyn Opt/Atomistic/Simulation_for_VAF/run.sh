#compute a short trajectory saving each time step velocity

gmx grompp -f ../mdp/prd_vaf.mdp -o prd_vaf -pp prd_vaf -po prd_vaf -c ../eql2 -t ../eql2 -p ../topol.top
gmx mdrun -deffnm prd_vaf -v

#compute VAF
gmx velacc -f prd_vaf.trr -s prd_vaf.tpr -n ../index.ndx -o vaf_atom.xvg

