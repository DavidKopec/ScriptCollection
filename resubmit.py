import os
import sys

import subprocess # maybe helpfull here ?

# def is_crashed(job_name:str, files:list[str]) ->bool:
#     files = [x for x in files if job_name in x]

#     for file in files:
#         if
#     pass

def resubmit(command:list[str],win_count:int) -> None:
    path:str =""
    for i in range(win_count):
        path = f"win_{i}"
        if(os.getcwd().split("/")[-1]!=path):
            os.chdir(os.getcwd()+"/"+path)
        print(os.getcwd())
        out = str(subprocess.check_output("pinfo"))

        if("Job is in error state!" in out):
            os.system("pkill -y")
            os.system("premovertf -y")
            os.system("psubmit "+" ".join(command))
        
        os.chdir(os.getcwd()+"/../")

    pass

def main() -> None:
    if(len(sys.argv)<3):
        print("Usage: win_cout cpu job_name ncpus=x walltime=1d -y")
        sys.exit()

    resubmit(sys.argv[2:],int(sys.argv[1]))

    

if __name__=="__main__":
    main()