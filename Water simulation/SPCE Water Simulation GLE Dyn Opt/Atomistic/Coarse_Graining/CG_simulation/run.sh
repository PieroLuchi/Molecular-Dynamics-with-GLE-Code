cp /home/piero/Scrivania/SPCE Water Simulation GLE Dyn Opt/Atomistic/Coarse_Graining/Force-maching/CG_CG_table.xvg /home/piero/Scrivania/SPCE Water Simulation GLE Dyn Opt/Atomistic/Coarse_Graining/CG_simulation/CG_CG_table.xvg

cp /home/piero/Scrivania/SPCE Water Simulation GLE Dyn Opt/Atomistic/Coarse_Graining/conf.gro /home/piero/Scrivania/SPCE Water Simulation GLE Dyn Opt/Atomistic/Coarse_Graining/CG_simulation/conf.gro

gmx grompp -f File/grompp.mdp -c conf.gro -n index.ndx -p topol.top 
gmx mdrun -table CG_CG_table.xvg
