#!/usr/bin/python
#MEANT TO MIMIC THE MCELL RUN OF THE OREGONATOR

import math
import os
from input import *

#in molecules
X = x_count
Y = y_count
Z = z_count

#rates in 1/s or 1/(molec*s)
ay_rate = k1*A  
xy_rate = k2 /LNa 
ax_rate = k3*A  
xx_rate = k4 / LNa
z_rate = k5*B  

# A + Y -> X             [ay_rate]
# X + Y -> C             [xy_rate]
# A + X -> X + X + 2Z    [ax_rate]
# X + X -> D             [xx_rate]
# Z + B -> fY            [z_rate]

xout = open("./output/xout",'w')
yout = open("./output/yout",'w')
xyout = open("./output/xyout",'w')
zyout = open("./output/zyout",'w')
zout = open("./output/zout",'w')


print"A2 = "+str(A)
print"B2 = "+str(B)
print"ay_rate = "+str(ay_rate)
print"ax_rate = "+str(ax_rate)
print"z_rate = "+str(z_rate)
print "step:\tX:\tY\tdx\tdy"

i=0

dx = 0
dy = 0
dz = 0

while i<iterations:
	i=i+1

	dx = (ay_rate*Y \
		-xy_rate*X*Y \
		+ax_rate*X \
		-2*xx_rate*X*X)*timestep
	dy = (-ay_rate*Y \
		-xy_rate*X*Y \
		+0.5*z_rate*Z)*timestep

	dz = (2*ax_rate*X \
		-z_rate*Z)*timestep

	X = X + dx
	Y = Y + dy
	Z = Z + dz

	xout.write(str(timestep*i*1e6)+'\t'+str(X) +'\n')
	yout.write(str(timestep*i*1e6)+'\t'+str(Y) +'\n')
	xyout.write(str(X)+'\t'+str(Y) +'\n')
	zout.write(str(i*timestep*1e6)+'\t'+str(Z) +'\n')
	zyout.write(str(Z)+'\t'+str(Y) +'\n')

	if i%10000000==0:
		print "step " + str(i) + "of " + str(iterations)


xout.close()
yout.close()
zout.close()
xyout.close()
zyout.close()

print "final count of X: " + str(X)
print "final count of Y: " + str(Y)
print "final count of Z: " + str(Z)




