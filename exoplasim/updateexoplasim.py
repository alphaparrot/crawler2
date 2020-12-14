import os

if __name__=="__main__":
    
    os.system('python -c "import exoplasim as exo; print(exo.__path__)" | tail -c +3 | head -c -3 > tmppath')
    with open('tmppath',"r") as spf:
        sourcedir = spf.read().strip()
        print("Found ExoPlaSim in %s"%sourcedir)
    os.system("rm tmppath")
    os.system("python -m pip install --user --upgrade --upgrade-strategy only-if-needed exoplasim")
    with open("%s/__init__.py"%sourcedir,"r") as sourcef:
        sourcecode = sourcef.read().split('\n')
    sourcecode[2] = 'sourcedir = "%s"'%sourcedir
    sourcecode = '\n'.join(sourcecode)
    #os.system("cp %s/__init__.py %s/preinit.py"%(sourcedir,sourcedir))
    try:
        with open("%s/__init__.py"%sourcedir,"w") as sourcef:
            sourcef.write(sourcecode)
        cwd = os.getcwd()
        os.chdir(sourcedir)
        os.system("source ~/.bashrc && prep_plasim && module load netcdf && module list "+
                  "&& ./configure.sh && cd postprocessor && ./build_init.sh")
        #os.system("./configure.sh")
        #os.system("cd postprocessor && ./build_init.sh")
        os.chdir(cwd)
    except PermissionError:
        raise PermissionError("\nHi! Welcome to ExoPlaSim. It looks like this is the first "+
                                "time you're using this program since installing, and you "+
                                "may have installed it to a location that needs root "+
                                "privileges to modify. This is not ideal! If you want to "+
                                "use the program this way, you will need to run python code"+
                                " that uses ExoPlaSim with sudo privileges; i.e. sudo "+
                                "python3 myscript.py. If you did this because pip install "+
                                "breaks without sudo privileges, then try using \n\n\tpip "+ "install --user exoplasim \n\ninstead. It is generally a "+
                                "very bad idea to install things with sudo pip install.")
