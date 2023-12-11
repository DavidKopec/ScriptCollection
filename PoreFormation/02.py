import os
import sys

def main(*args):
    
    file_count:int = 0
    path:str=""    
    
    if(len(sys.argv)>1):
        file_count=int(sys.argv[1])
        #input(f"fc: {file_count}")
    
    # with open("tpr_files.dat","x") as f:
    #     for i in range(1,file_count+1):
    #         path = f"win_{i}"
    #         f.write(f"{path}/storage/md0001.tpr\n")
    
    fill:str = ""
    KAPPA:float
    AT:float

    with open("metadata.dat","x") as f:
        for i in range(file_count+1):
            fill=""
            path = f"win_{i}"
            
            for j in range(1,5):
                with open(f"{path}/storage/md000{j}.win-55.dat","r") as px:
                    file = px.readlines()[1:]
                    for line in file:
                        fill+=f"""{line.split(" ")[1]} {line.split(" ")[2]}\n"""
            
            with open(f"{path}/storage/pullx.xvg","x") as c:
                c.write(fill)

            with open(f"{path}/plumed.dat","r") as pl:
                content = pl.readlines()

                for line in content:
                    if("us_tails:" in line):
                        KAPPA = float(line.split()[3].split('=')[1])
                        AT = float(line.split()[4].split('=')[1])
                
                fill = "%f\t %f\n" % (AT,KAPPA)


            f.write(f"{path}/storage/pullx.xvg\t{fill}")

    #os.system("./wham ")


if __name__=="__main__":
    main()
