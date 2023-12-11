import sys
import os

def win_reduction(start:int,stop:int,step:int = 0):
    wins = [file for file in os.listdir() if "win-" in file]
    
    if step==0: step=len(wins)

    for x,win in enumerate(wins):
        if(start<=x<=stop):
            if(x%step==0):
                os.remove(win)
                continue
            os.remove(win)


def main():
    args = " ".join(sys.argv).split()

    if(len(args)>2):
        if(len(args)==4):
            win_reduction(args[1],args[2],args[3])
            return
        win_reduction(args[1],args[2])
    else:
        print("Usage: python3 ReduceWindows.py Start Stop [Step]")

if __name__=="__main__":
    main()