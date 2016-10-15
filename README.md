# ptimsa
add molecules to areas with high number density for molecular dynamics configurations


tested with gfortran 5.4.0 on cygwin and 4.4.7 on linux


compile with 'gfortran -O3 ptimsa.f95 -o ptimsa.e'


example usage:


./ptimsa -c spc216.gro -ci phenol.pdb -n 5 -g SOL -dl SOL -oc phenol-in-water.pdb 
