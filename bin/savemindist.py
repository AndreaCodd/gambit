import numpy as np
from esys.escript import length
#
csvname = "mindist.csv"

# import measuring locations
pts=np.loadtxt("MtIsaCloncurry_MeasEs.csv",delimiter=',')
mdist=[]
for ind1 in range(len(pts)):
    dist = 5000 
    for ind2 in range(len(pts)):
        if ind1 != ind2:
            dist2= length(pts[ind1]-pts[ind2])
            if dist2 < dist :
                dist = dist2
    mdist.append(dist)
mdist=np.array(mdist)
np.savetxt("MtIsaCloncurry_mindist.csv",mdist, delimiter=",")
