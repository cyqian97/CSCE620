from DelVorCov import *
from time import time
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# Random vertices
n=15
seed = int(time())
np.random.seed(seed)
print("Random seed: {}".format(seed))
vcoords = np.random.uniform(0,1,(n,3))

# # n=5 diamond
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

# # n = 12 regular icosahedron
# phi = (5**0.5-1)/2
# vcoords = np.array([
#     [0,1,phi],
#     [0,1,-phi],
#     [0,-1,phi],
#     [0,-1,-phi],
#     [1,phi,0],
#     [1,-phi,0],
#     [-1,phi,0],
#     [-1,-phi,0],
#     [phi,0,1],
#     [phi,0,-1],
#     [-phi,0,1],
#     [-phi,0,-1],
# ])

start_time = time()
db = DelVorCov(vcoords)
print("--- %s seconds ---" % (time() - start_time))
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

ax = plt.figure().add_subplot(projection='3d')
verts = []
for bid in db.convhull.all_bids:
    verts.append(db.getBoundaryVCoords(bid))

poly = Poly3DCollection(verts, linewidths=1, alpha=0.2)
poly.set_facecolor("tab:blue")
poly.set_edgecolor("black")
ax.add_collection3d(poly)
ax.scatter(vcoords[:, 0], vcoords[:, 1], vcoords[:, 2])

# vid = 4
# p1 = db.vcoords[vid]
# for eid in db.v2e[vid]:
#     if db.e2v[eid][0] == vid:
#         p2 = db.vcoords[db.e2v[eid][1]]
#     else:
#         p2 = db.vcoords[db.e2v[eid][0]]
#     ax.plot3D([p1[0],p2[0]], [p1[1],p2[1]],[p1[2],p2[2]], 'tab:red')

ax.quiver(db.vcoords[db.convhull.all_vids,0],
          db.vcoords[db.convhull.all_vids,1],
          db.vcoords[db.convhull.all_vids,2],
          db.convhull.v_norms[db.convhull.all_vids,0],
          db.convhull.v_norms[db.convhull.all_vids,1],
          db.convhull.v_norms[db.convhull.all_vids,2],
          length = 0.5)
ax.axis('equal')   
plt.show()