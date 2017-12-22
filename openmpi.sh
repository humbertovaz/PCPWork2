# load do modulo


module load gnu/openmpi_eth/2.0.0
cd /home/a73236/PCPWork2
export FILE=heatplatempi.c
mpicc $FILE -o teste
mpirun -np $1 --mca btl tcp,sm,self ./teste $2 $3 $4 $5 $6


