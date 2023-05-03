from Geometry import *
from time import time
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

class DelVorCov(GeometryDatabase):
    def __init__(self, vcoords):
        super().__init__(vcoords)
        self.boywerWatson()

    def boywerWatson(self):
        d = 100
        far_points = np.array([[0, 0, d], [d, 0, -d], [-d, d, -d], [-d, -d, -d]])
        vcoords = self.vcoords
        far_ids = [len(vcoords),len(vcoords)+1,len(vcoords)+2,len(vcoords)+3]
        vcoords = np.block([[vcoords],[far_points]])
        db = GeometryDatabase(vcoords)
        current_simplices = [Simplex(far_ids, db)]
        for vid in range(len(vcoords)-4):
            v = vcoords[vid]
            bad_simplices = []
            pop_ids = []
            for sid, s in enumerate(current_simplices):
                if s.isInCircumcicle(v):
                    bad_simplices.append(s)
                    pop_ids.append(sid)
            pop_ids.reverse()
            for sid in pop_ids:
                current_simplices.pop(sid)

            bad_boundaries = []
            for s in bad_simplices:
                bad_boundaries += s.bids

            # Counting sort
            bid_max = np.max(bad_boundaries)
            bid_count = np.zeros(bid_max+1)
            for bid in bad_boundaries:
                bid_count[bid] += 1

            for bid in range(bid_max+1):
                if bid_count[bid] == 1:
                    current_simplices.append(
                        Simplex(list(db.b2v[bid]) + [vid], db))

        pop_ids = []
        for sid, s in enumerate(current_simplices):
            if len(set.intersection(set(far_ids), s.vids)) > 0:
                pop_ids.append(sid)
        pop_ids.reverse()
        for sid in pop_ids:
            current_simplices.pop(sid)
        
        self.v2b = db.v2b[:-4]
        self.b2v = db.b2v
        self.delaunay = current_simplices

        return

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
current_simplices = db.delaunay
print("Num Simplices: {}".format(len(current_simplices)))
all_boundaries = []
for s in current_simplices:
    all_boundaries += s.bids

ax = plt.figure().add_subplot(projection='3d')
verts = []
for b in set(all_boundaries):
    vids = db.b2v[b]
    verts.append(db.getBoundaryVCoords(b))

poly = Poly3DCollection(verts, linewidths=1, alpha=0.1)
poly.set_facecolor("tab:blue")
poly.set_edgecolor("black")
ax.add_collection3d(poly)
ax.scatter(vcoords[:, 0], vcoords[:, 1], vcoords[:, 2])
plt.show()