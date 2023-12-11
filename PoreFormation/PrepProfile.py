import os 

coor=0.0
free=0.0
fill=""

with open("G.xvg","r") as G:
	for line in G.readlines()[1:]:
		if('#' in line):
			break
		coor = float(line.split()[0])
		free = float(line.split()[1])
		
		fill+= "%f %f\n" % (coor,free)

with open("G_rdy.xvg","w") as G:
	G.write(fill)


		
