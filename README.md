# ptimsa
add molecules to areas with high number density for molecular dynamics configurations


tested with gfortran 5.4.0 on cygwin and 4.4.7 on linux


compile with 'gfortran -O3 catnip.f95 -o catnip.e'


example usage:


./catnip.e -c spc216.gro -ci phenol.pdb -n 5 -g SOL -dl SOL -oc phenol-in-water.pdb 


running ./catnip.e with -h or nothing will show more options
