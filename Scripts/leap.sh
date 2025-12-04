#!/bin/bash

#SBATCH -J 1-leap
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --partition=gpu
#SBATCH --mem=10G

export OMP_NUM_THREADS=$SLURM_NPROCS
module load amber/24

for i in sistem
do

cd $i

echo -e "
source leaprc.protein.ff19SB
source leaprc.gaff
source leaprc.water.opc

set default PBRadii mbondi2

com = loadpdb $i.pdb

check com
charge com
addions com Na+ 40
addions com Cl- 30
solvateBox com OPCBOX 10.0
saveamberparm com $i.prmtop $i.inpcrd
savepdb com Sistema-$i.pdb
savemol2 com Sistema-$i.mol2 1
quit

" > leap.in

tleap -f leap.in

cd ..

done
