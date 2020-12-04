import glob
import os
import sys

template = ("#!/bin/bash -l   \n"+
		"#PBS -l nodes=1:ppn=1  \n"+
		"#PBS -q workq  \n"+
		"#PBS -m a     \n"+
		"#PBS -r n     \n"+
		"#PBS -l walltime=03:00:00        \n"+
		"#PBS -N transit_%s                   \n"+
		"# EVERYTHING ABOVE THIS REQ'D    \n"+
		"cd $PBS_O_WORKDIR       \n"+
		"module load gcc/4.9.1   \n"+
		"module load python/2.7.9 \n"+
		"mkdir /mnt/node_scratch/paradise/%s     \n"+
		"cp %s /mnt/node_scratch/paradise/%s/    \n"+
                "cp fulltransit.py /mnt/node_scratch/paradise/%s/   \n"+
		"cd /mnt/node_scratch/paradise/%s/   \n"+
		"python fulltransit.py %s   \n"+
                "cp *.png $PBS_O_WORKDIR/   \n"+
                "cp *.pdf $PBS_O_WORKDIR/   \n"+
                "cp *.npy $PBS_O_WORKDIR/   \n"+
                "rm -rf *   \n"+
                "cd $PBS_O_WORKDIR   \n")


if __name__=="__main__":
  if len(sys.argv)==1:
  	files = glob.glob("*.nc")
  	for f in files:
     		with open("run"+f[:-3],"w") as runf:
        		runf.write(template%(f[:-3],f[:-3],f,f[:-3],f[:-3],f[:-3],f))
     		os.system("qsub run"+f[:-3])
  else:
	for f in sys.argv[1:]:
		with open("run"+f[:-3],"w") as runf:
                        runf.write(template%(f[:-3],f[:-3],f,f[:-3],f[:-3],f[:-3],f))
                os.system("qsub run"+f[:-3])
