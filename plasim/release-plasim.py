import os
import sys
import glob

if __name__=="__main__":
    outputfiles = sorted(glob.glob("*.nc"))
    if "all" not in sys.argv[:]:
        newest = -1
        for f in outputfiles:
            age = os.path.getmtime(f)
            if age>newest and os.path.getsize(f)>1.0e6:
                newest=age
                lastfile = f
        name = sys.argv[1]
        os.system("cp "+lastfile+" ../output/"+name+".nc")
    else:
        name = sys.argv[1]
        os.system("mkdir ../output/"+name)
        for f in outputfiles:
            os.system("cp "+f+" ../output/"+name+"/"+name+"_"+f)
        os.system("cp *DIAG* ../output/"+name+"/")
    os.system("rm -rf *.nc")
    os.system("rm -rf MOST*")
    
    os.system("python release.py")