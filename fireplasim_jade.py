import numpy as np

if __name__=="__main__":
        header = "# PID MODEL JOBNAME STATUS NCORES QUEUE nfixorb eccen obliq vernlon lockedyear year naqua pCO2 pressure flux script extra notify alloutput restart"

	tstr1 = " 0 8 greenq 1 0.0 0.0 90.0 10.0/0.0 36 1 363.96 " 
	tstr2 = "equilibriate.sh relax.py a 0 locked_template"
	tstr3 = "equilibriate.sh relax.py a 0 locked_frozen"

	ttxtf = open("tasks.crwl","r")
	ttxt = ttxtf.read()
	ttxtf.close()
        ttxt += header + "\n"
	n=3103

        pc0 = 3.6396e-4

	for sol in np.linspace(700,1500,num=25):
		for ps in 10**np.linspace(np.log10(0.2),np.log10(20.0),num=50):
			ttxt+=str(n)+" plasim psgrid"+str(n)+tstr1+str(ps+pc0)+" "+str(sol)+" "+tstr2+'\n'
			n+=1
			ttxt+=str(n)+" plasim pscgrid"+str(n)+tstr1+str(ps+pc0)+" "+str(sol)+" "+tstr3+'\n'
			n+=1
	
	ttxtf = open("tasks.crwl","w")
	ttxtf.write(ttxt)
	ttxtf.close()