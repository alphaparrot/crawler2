import numpy as np
import glob
import os, sys
import exoplasim as exo


gplasim = True
TIMELIMIT = 1440.0*1.75

#This version lets the model relax.


if __name__=="__main__":
  if gplasim:
    if not os.path.exists("weathering.pso"):
        wf=open("weathering.pso","w")
        wf.write("     CO2       AVG SURF T   WEATHERING    OUTGASSING      DpCO2       NEW CO2\n")
        wf.close()

  os.system("rm keepgoing")
  
  os.system("rm -f Abort_Message")
  os.system("echo 'SURFACE      TOA'>balance.log")
  os.system("echo 'SURFACE      TOA'>slopes.log")
  
  try:
    model = np.load("model.npy",allow_pickle=True).item()
  except:
    model = np.load("model.npy").item()
  modelodir = model.workdir
  model.workdir = os.getcwd()
  
  equilibrium = model.runtobalance(threshold = model.threshold,
                                   timelimit = TIMELIMIT)

  model.workdir = modelodir
  model.save()
  if not equilibrium:
      os.system("touch keepgoing")