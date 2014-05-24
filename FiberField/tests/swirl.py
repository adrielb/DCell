import numpy  as np
import pandas as pd
import matplotlib
matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt
import os
from matplotlib.patches import Rectangle
import matplotlib.animation as animation

petscdir = os.environ[ 'PETSC_TMP' ]
os.chdir( petscdir )

'''
plt.ion()
'''

params = pd.read_csv( "params.csv", index_col=[0,1] )
print params

ranks = params.index.get_level_values( "rank" ).unique()

df = pd.read_csv("x.csv", index_col=[0,1])
df.head(19)

df.loc[0]


vIDs = df.index.get_level_values( 'vID' ).unique()
vIDs.sort()
vCol = np.random.rand(vIDs.size,3)

#bug in pandas getting max without unique
timax = df.index.get_level_values( 'ti' ).unique().max()

ti = 2
fig = plt.figure(1)
fig.clear()
ax = fig.gca()
ax.set_title("")
ax.set_xlabel("")
ax.set_ylabel("")
ax.set_xlim(-1, 1)
ax.set_ylim(-1, 1)
ax.plot(df['x'][0],df['z'][0], label=''  ,marker='o', lw=0)
#ax.plot(df['x'][1],df['z'][1], label=''  ,marker='o', lw=0)
#ax.plot(df['x'][2],df['z'][2], label=''  ,marker='o', lw=0)
#ax.legend(loc='upper left')
fig.canvas.draw()
plt.show()

def drawBBox( ax, rank, var, **kwargs ):
    minmax = params[["x","z"]].loc[rank]
    lmin = minmax.loc[var+".min"]
    w = minmax["x"].loc[var+".max"] - lmin["x"]  
    h = minmax["z"].loc[var+".max"] - lmin["z"]  
    ax.add_patch( Rectangle( lmin, w, h, **kwargs ) )

def drawPartitions( ax ):
    drawBBox( ax, 0, "globalBounds", lw=10, edgecolor="red", facecolor="none" )
    for r in ranks:
        drawBBox( ax, r, "localBounds" , edgecolor="black", facecolor="none" ) 

fig = plt.figure(2)
fig.clear()
ax = fig.gca()
ax.set_title("")
ax.set_xlabel("")
ax.set_ylabel("")
ax.set_xlim(-1, 1)
ax.set_ylim(-1, 1)
drawPartitions(ax)
for vID in vIDs:
    X = df.xs(vID, level='vID' ).sort()
    ax.plot(X["x"], X["z"], label=''  ,marker='o', lw=1)
#ax.legend(loc='upper left')
fig.canvas.draw()
fig.show()

fig=plt.figure(3)
fig.clear()
ax=fig.gca()
ax.set_xlim(-1, 1)
ax.set_ylim(-1, 1)
drawPartitions(ax)
pts = ax.scatter(df['x'][0],df['z'][0],s=40,c=vCol)
#pts = ax.plot([],[],lw=0,marker='.')
fig.canvas.draw()
fig.show()

def myinit():
    #pts[0].set_array([],[])
    #pts.set_array(np.array([]))
    pts.set_offsets(df[['x','z']].loc[0].values)
    return pts
def myani( i ):
    #pts[0].set_data(df['x'][i],df['z'][i])
    pts.set_offsets(df[['x','z']].loc[i].values)
    return pts

ani = animation.FuncAnimation(fig, myani, np.arange(0,timax),
        interval=200, blit=True, init_func=myinit )
fig.show() 

############3

fig=plt.figure(3)
ax=fig.gca()
ax.set_xlim(-1, 1)
ax.set_ylim(-1, 1)
for i in np.arange( 0,timax):
    fig.clear()
    drawPartitions(ax)
    pts = ax.scatter(df['x'][i],df['z'][i],s=40,c=vCol)
    plt.draw()
    plt.pause(1)
    plt.clf()
