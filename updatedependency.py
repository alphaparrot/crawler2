import os,sys

if __name__=="__main__":
   pid = int(sys.argv[1])
   with open("%d_new.id"%pid,"r") as nf:
       newid = nf.read().split('\n')[0].split()[0]
   with open("%d.id"%pid,"r") as oldf:
       children = oldf.read().split('\n')[1:]
   
   for child in children:
       os.system("qalter -W depend=afterany:%s %s"%(newid,child))
   
   os.system("mv %d_new.id %d.id"%(pid,pid))
   

