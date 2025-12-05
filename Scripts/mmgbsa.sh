PRO=_prod
System=S_6

echo "
&general
startframe=1500,
endframe=2000,
interval=1,
receptor_mask=":34-198",
ligand_mask=":1-33",
strip_mask="!:1-198",
keep_files=0
/
&gb
   igb=5
   saltcon=0.10
/
" > mmgbsa.in

ante-MMPBSA.py -p $System.prmtop -s "!:1-198" -c complex.prmtop --radii=mbondi2
ante-MMPBSA.py -p complex.prmtop -n ":1-33" -r receptor.prmtop -l ligand.prmtop --radii=mbondi2
MMPBSA.py -O -i mmgbsa.in -o MMGBSA_M4_$System.dat -sp $System.prmtop -cp complex.prmtop -rp receptor.prmtop -lp ligand.prmtop -y $System$PRO.nc
 
