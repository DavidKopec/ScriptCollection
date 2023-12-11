import os
import sys

import subprocess # maybe helpfull here ? YES !!!!

def resubmit(command:list[str],win_count:int) -> None:
    path:str =""
    for i in range(win_count):
        path = f"win-{win_count}"
        os.chdir(os.getcwd()+"/"+path)
        out = subprocess.check_output("pinfo")

        if("Job is in error state!" in out):
            os.system("psubmit"+" ".join(command))

    pass

def main() -> None:
    if(len(sys.argv)<3):
        print("Usage: win_cout cpu job_name ncpus=x walltime=1d -y")
        sys.exit()

    resubmit(sys.argv[2:],sys.argv[1])

    

if __name__=="__main__":
    main()