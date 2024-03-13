#!/bin/bash
#SBATCH --time=2:00:00                     # set max wallclock time
#SBATCH --output=job-%N-%j.out             # set output log name
#SBATCH --error=job-%N-%j.error            # set error log name %N=first node, %j=job number
#SBATCH --job-name=isolated-pot	           # set job name
#SBATCH -p all                             # set partition (queue)
#SBATCH --nodes=1                          # set node number
#SBATCH --ntasks-per-node=16               # set task number
#SBATCH --mem-per-cpu=2G
#SBATCH --mail-type=ALL                    # set mail alerts
#SBATCH --mail-user=simone.gordon@ed.ac.uk # set mail address

. /usr/share/modules/init/bash

cd ${SLURM_SUBMIT_DIR}

module load hdf5/gcc-11/mpich/1.14.0
export LD_LIBRARY_PATH=/home/sgordon/software/local/lib:$LD_LIBRARY_PATH

#./simrun.pl -mpi "mpirun -np 18"
mpirun -np 16 ./enzo.exe -d NohProblem3D.enzo >& estd.out 2>&1
