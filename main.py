from DelVorCov import *
from time import time
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# Random vertices
# n=20
# seed = int(time())
# np.random.seed(seed)
# print("Random seed: {}".format(seed))
# vcoords = np.random.uniform(0,1,(n,3))

# n=5 diamond
# vcoords = np.array([[0,0,1],
#           [np.cos(np.pi/6),np.sin(np.pi/6),0],
#           [np.cos(5*np.pi/6),np.sin(5*np.pi/6),0],
#           [0,-1,0],
#           [0,0,-1]])*1.5

# n=6 prism
# vcoords = np.array([[np.cos(np.pi/6), np.sin(np.pi/6), 1],
#                     [np.cos(5*np.pi/6), np.sin(5*np.pi/6), 1],
#                     [0, -1, 1],
#                     [np.cos(np.pi/6), np.sin(np.pi/6), -1],
#                     [np.cos(5*np.pi/6), np.sin(5*np.pi/6), -1],
#                     [0, -1, -1]])

# n = 12 regular icosahedron
phi = (5**0.5-1)/2
vcoords = np.array([
    [0,1,phi],
    [0,1,-phi],
    [0,-1,phi],
    [0,-1,-phi],
    [1,phi,0],
    [1,-phi,0],
    [-1,phi,0],
    [-1,-phi,0],
    [phi,0,1],
    [phi,0,-1],
    [-phi,0,1],
    [-phi,0,-1],
])


db = DelVorCov(vcoords)
# Visualization

ax = plt.figure().add_subplot(projection='3d')
verts = []
for bid in db.all_bids:
    verts.append(db.getBoundaryVCoords(bid))

poly = Poly3DCollection(verts, linewidths=1, alpha=0.1)
poly.set_facecolor("tab:blue")
poly.set_edgecolor("black")
ax.add_collection3d(poly)
ax.scatter(vcoords[:, 0], vcoords[:, 1], vcoords[:, 2])
plt.show()