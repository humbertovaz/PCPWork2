# load do modulo


module load gnu/openmpi_eth/2.0.0
cd /home/a73236/PCPWork2
export FILE=heatplatempi.c
mpicc $FILE -o teste
./teste


