import os
import sys

winCount:int = 0
skip:int = -1

if(len(sys.argv)>1):
    winCount = int(sys.argv[1])
    skip = int(sys.argv[2])

mainDir = os.getcwd()

for i in range(0,winCount+1):
    
    if(i==skip):
        continue

    winName = ""
    
    winName = "win_"+str(i)
    os.system(f"cp precycle_md_FE ./{winName}")
    
    os.chdir(mainDir+"/"+winName)
    #os.system("mv membrane.top system.top")
    os.system("premovertf -f")
    os.system("psubmit -y cpu precycle_md_FE ncpus=8, walltime=1d")
