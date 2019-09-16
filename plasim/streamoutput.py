import numpy as np
import sys

def readchunk(f,nrows,ncols):
    data = np.zeros((nrows,ncols))
    for n in range(nrows):
        try:
            line = f.next()
            line = line.split()
            for k in range(ncols):
                data[n,k] = float(line[k])
        except:
            return data[:n,:],False
    return data,True


if __name__=="__main__":
    
    fname = sys.argv[1]
    ncols = int(sys.argv[2])
    try:
        baseline = int(sys.argv[3])
    except:
        baseline = 23040 #One month with 32 timesteps per day
    
    output = np.zeros((ncols,1))
    
    with open(fname,"r") as fstream:
        keepgoing = True
        while keepgoing:
            chunk,keepgoing = readchunk(fstream,baseline,ncols)
            chunk2 = np.zeros((ncols,1))
            chunk2[:,0] = np.nanmean(chunk,axis=0)
            output = np.append(output,chunk2,axis=1)
        
    output = output[:,1:]
    
    np.save(fname[:-4]+".npy",output)
    
    
    