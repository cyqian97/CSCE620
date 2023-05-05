from DelVorCov import *
from GenPointCloud import *
from time import time
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

n = 20
vcoords = genPointCloud("rand",n=n)
ax = plt.figure().add_subplot(projection='3d')
ax.scatter(vcoords[:, 0], vcoords[:, 1], vcoords[:, 2],s = 2)
 
ax.view_init(elev=25, azim=-40, roll=0) 
ax.set_xlim(-2,1)
ax.set_ylim(-1.5,1.5)
ax.set_zlim(-1,2)
plt.show()
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

ax.quiver(db.vcoords[db.convhull.all_vids,0],
          db.vcoords[db.convhull.all_vids,1],
          db.vcoords[db.convhull.all_vids,2],
          db.convhull.v_norms[db.convhull.all_vids,0],
          db.convhull.v_norms[db.convhull.all_vids,1],
          db.convhull.v_norms[db.convhull.all_vids,2],
          length = 0.5)
ax.axis('equal')   
plt.show()

start_time = time()
db2, crust = db.Amenta()
print("--- %s seconds ---" % (time() - start_time))

print(len(crust))

# ax = plt.figure().add_subplot(projection='3d')
# verts = []
# for bid in db2.all_bids:
#     vert = db2.getBoundaryVCoords(bid)
#     if np.min(vert[:,2])>=0:
#         verts.append(vert)
#     # verts.append(db2.getBoundaryVCoords(bid))

# poly = Poly3DCollection(verts, linewidths=1, alpha=0.1)
# poly.set_facecolor("tab:blue")
# poly.set_edgecolor("black")
# ax.add_collection3d(poly)
# ax.scatter(vcoords[:, 0], vcoords[:, 1], vcoords[:, 2])
# plt.show()


ax = plt.figure().add_subplot(projection='3d')
verts = []
for bid in crust:
    vert = db2.getBoundaryVCoords(bid)
    verts.append(vert)
    # if np.min(vert[:,2])>0:
    #     verts.append(vert)

poly = Poly3DCollection(verts, linewidths=1, alpha=0.2)
poly.set_facecolor("tab:blue")
poly.set_edgecolor("black")
ax.add_collection3d(poly)
ax.scatter(vcoords[:, 0], vcoords[:, 1], vcoords[:, 2],s = 0.2)

ax.view_init(elev=25, azim=-40, roll=0) 
ax.set_xlim(-2,1)
ax.set_ylim(-1.5,1.5)
ax.set_zlim(-1,2)

plt.show()

# all_correct = []
# for i in np.sqrt(n):
#     for j in np.sqrt(n):

