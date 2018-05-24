import numpy as np

if __name__=="__main__":
        header = "# PID MODEL JOBNAME STATUS NCORES QUEUE nfixorb eccen obliq vernlon lockedyear year naqua flux gravity radius pCO2 pressure script extra notify alloutput restart"

	tstr1 = " 0 8 sandyq 1 0.0 0.0 90.0 11.2/0.0 36 1 956.0 10.9 1.1 " 
	tstr2 = "equilibriate.sh relax.py a 0 locked_restart"

	ttxtf = open("tasks.crwl","r")
	ttxt = ttxtf.read()
	ttxtf.close()
        ttxt += header + "\n"
	n=1

	for pc in 10**np.linspace(-2,6,num=60):
		for ps in 10**np.linspace(-1,np.log10(20.0),num=50):
			ttxt+=str(n)+" plasim pgrid"+str(n)+tstr1+str(pc)+" "+str(ps+pc*1.0e-6)+" "+tstr2+'\n'
			n+=1
	
	ttxtf = open("tasks.crwl","w")
	ttxtf.write(ttxt)
	ttxtf.close()