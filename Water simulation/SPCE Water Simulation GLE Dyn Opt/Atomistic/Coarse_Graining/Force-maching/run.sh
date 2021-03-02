#! /bin/bash -e

#equilibration time in Gromacs units (ps)
equi=1
echo equi = $equi

if [[ ! -f ../../prd.tpr || ! -f ../../prd.trr ]]; then
  echo "Run atomistic simulation first"
  exit 1
fi


echo "Running force matching"
csg_fmatch --top ../../prd.tpr --trj ../../prd.trr --begin $equi  --options fmatch.xml --cg water.xml --nframes 500 --verbose

#integrate force table to get potential
csg_call table integrate CG-CG.force CG-CG.pot
csg_call table linearop CG-CG.pot CG-CG.pot -1 0


#copy CG-CG.pot to new file to prevent deletion by clean command
cp CG-CG.pot input.pot

#convert to gromacs potential
csg_call --options fmatch.xml --ia-name CG-CG --ia-type non-bonded convert_potential gromacs --clean input.pot table_CG_CG.xvg

