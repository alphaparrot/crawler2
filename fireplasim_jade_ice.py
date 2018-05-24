import numpy as np

if __name__=="__main__":
        header = "# PID MODEL JOBNAME STATUS NCORES QUEUE nfixorb eccen obliq vernlon lockedyear year naqua pCO2 pressure flux script extra notify alloutput keeprestart restart"

	tstr1 = " 0 8 hpq 1 0.0 0.0 90.0 10.0/0.0 36 1 363.96 " 
	tstr2 = "super-equil.sh super_relax.py a 1 1 "
	tstr3 = "super-equil.sh super_relax.py a 1 1 "

	ttxtf = open("tasks.crwl","r")
	ttxt = ttxtf.read()
	ttxtf.close()
        ttxt += header + "\n"
	n=10
	l=33
	k=151

        ps = 1.011
        sol = 1050.0
        for p in range(0,40):
            ttxt+=(str(k)+" plasim jade_ice"+str(k)+tstr1+str(ps)+" "+str(sol)
                    +" "+tstr2+'jade'+str(n)+'_restart \n')
            n+=1
            k+=1
        for p in range(0,18):
            ttxt+=(str(k)+" plasim jade_ice"+str(k)+tstr1+str(ps)+" "+str(sol)
                    +" "+tstr3+'cjade'+str(l)+'_restart \n')
            l+=1
            k+=1
            
        
        n=10
        l=33
        
        sol = 1200.0
        for p in range(0,40):
            ttxt+=(str(k)+" plasim jade_ice"+str(k)+tstr1+str(ps)+" "+str(sol)
                    +" "+tstr2+'jade'+str(n)+'_restart \n')
            n+=1
            k+=1
        for p in range(0,18):
            ttxt+=(str(k)+" plasim jade_ice"+str(k)+tstr1+str(ps)+" "+str(sol)
                    +" "+tstr3+'cjade'+str(l)+'_restart \n')
            l+=1
            k+=1
            
	ttxtf = open("tasks.crwl","w")
	ttxtf.write(ttxt)
	ttxtf.close()