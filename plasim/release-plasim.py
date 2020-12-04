import os
import sys
import glob

if __name__=="__main__":
    outputfiles = sorted(glob.glob("*.nc"))
    os.system("mkdir ../output")
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
        newest = -1
        for f in outputfiles:
            age = os.path.getmtime(f)
            if age>newest and os.path.getsize(f)>1.0e6:
                newest=age
                lastfile = f
        name = sys.argv[1]
        os.system("cp "+lastfile+" ../output/"+name+".nc")
        
    restartfiles = sorted(glob.glob("*REST.*"))
    if "restart" in sys.argv[:]:
        if "all" in sys.argv[:]:
            for f in restartfiles:
                os.system("cp "+f+" ../output/"+name+"/"+name+"_"+f)
        
    snapshotfiles = sorted(glob.glob("snapshots/*.nc"))
    if "all" not in sys.argv[:]:
        newest = -1
        for f in snapshotfiles:
            age = os.path.getmtime(f)
            if age>newest and os.path.getsize(f)>1.0e6:
                newest=age
                lastfile = f
        name = sys.argv[1]+"_snapshot"
        os.system("cp "+lastfile+" ../output/"+name+".nc")
    else:
        name = sys.argv[1]
        for s in snapshotfiles:
            f = s.split('/')[1]
            os.system("cp "+s+" ../output/"+name+"/"+name+"_"+f)
        newest = -1
        for f in snapshotfiles:
            age = os.path.getmtime(f)
            if age>newest and os.path.getsize(f)>1.0e6:
                newest=age
                lastfile = f
        name = sys.argv[1]+"_snapshot"
        os.system("cp "+lastfile+" ../output/"+name+".nc")
        
    hcfiles = sorted(glob.glob("highcadence/*.nc"))
    if "all" not in sys.argv[:]:
        newest = -1
        for f in hcfiles:
            age = os.path.getmtime(f)
            if age>newest and os.path.getsize(f)>1.0e6:
                newest=age
                lastfile = f
        name = sys.argv[1]+"_highcadence"
        if len(hcfiles)>0:
            os.system("cp "+lastfile+" ../output/"+name+".nc")
    else:
        name = sys.argv[1]
        for s in hcfiles:
            f = s.split('/')[1]
            os.system("cp "+s+" ../output/"+name+"/"+name+"_"+f)
        newest = -1
        for f in hcfiles:
            age = os.path.getmtime(f)
            if age>newest and os.path.getsize(f)>1.0e6:
                newest=age
                lastfile = f
        name = sys.argv[1]+"_highcadence"
        if len(hcfiles)>0:
            os.system("cp "+lastfile+" ../output/"+name+".nc")
        
    psofiles = sorted(glob.glob("*.pso"))
    if "all" not in sys.argv[:]:
        newest = -1
        lastfile = ''
        for f in psofiles:
            age = os.path.getmtime(f)
            if age>newest:
                newest=age
                lastfile = f
        name = sys.argv[1]
        os.system("cp "+lastfile+" ../output/"+name+".pso")
        os.system("cp cyclelog.txt ../output/"+name+"_cyclelog.txt")
    else:
        name = sys.argv[1]
        os.system("mkdir ../output/"+name)
        for f in psofiles:
            os.system("cp "+f+" ../output/"+name+"/"+name+"_"+f)
        os.system("cp *DIAG* ../output/"+name+"/")
        os.system("cp cyclelog.txt ../output/"+name+"/"+name+"_cyclelog.txt")
        
    hifiles = sorted(glob.glob("*.STORM"))
    if len(hifiles)>0:
        if "all" not in sys.argv[:]:
            newest = -1
            lastfile = ''
            for f in hifiles:
                age = os.path.getmtime(f)
                if age>newest:
                    newest=age
                    lastfile = f
            name = sys.argv[1]
            os.system("cp "+lastfile+" ../output/"+name+".STORM")
        else:
            name = sys.argv[1]
            os.system("mkdir ../output/"+name)
            for f in hifiles:
                os.system("cp "+f+" ../output/"+name+"/"+name+"_"+f)
            newest = -1
            lastfile = ''
            for f in hifiles:
                age = os.path.getmtime(f)
                if age>newest:
                    newest=age
                    lastfile = f
            name = sys.argv[1]
            os.system("cp "+lastfile+" ../output/"+name+".STORM")
                
        
    os.system("rm -rf *.pso")
    os.system("rm -rf *.nc")
    os.system("rm -rf snapshots/")
    os.system("rm -rf highcadence/")
    os.system("rm -rf MOST*")
    os.system("rm -rf cyclelog.txt")
    os.system("cp ../../release.py .")
    os.system("python release.py")