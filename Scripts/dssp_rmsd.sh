System=S_6
PRO=_prod

# Create input for cpptraj
cat > traj.in << EOF
trajin ${System}${PRO}.nc
strip :WAT,Na+,Cl-
autoimage origin :166-198
rms fit :166-198

reference ${System}.inpcrd

# RMSD per monomer (6 monomers)
rms M4  :1-33@CA      reference out rmsd_M1.dat mass
rms M5  :34-66@CA     reference out rmsd_M2.dat mass
rms M6  :67-99@CA     reference out rmsd_M3.dat mass
rms M7  :100-132@CA   reference out rmsd_M4.dat mass
rms M8  :133-165@CA   reference out rmsd_M5.dat mass
rms M9  :166-198@CA   reference out rmsd_M6.dat mass

# Fibril secondary structure
secstruct :1-198   out secstruct_fibril.dat \
                   sumout ss_fibril_sum.dat \
                   assignout ss_fibril_assign.dat

# Peptide secondary structure
secstruct :196-226 out secstruct_peptide.dat \
                   sumout ss_peptide_sum.dat \
                   assignout ss_peptide_assign.dat
EOF

# Run cpptraj
cpptraj ${System}.prmtop traj.in

print(f"Missing replicas for {system}")


