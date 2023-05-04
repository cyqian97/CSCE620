from Geometry import *


class DelVorCov(GeometryDatabase):
    def __init__(self, vcoords):
        super().__init__(vcoords)
        self.delaunay = [] # list of simplices
        # self.b2s = [] # boundaries to simplices
        self.convhull = self.ConvexHull() # list of convex hull bids
        self.boywerWatson()
        self.processConvHull()

    class ConvexHull():
        def __init__(self):
            self.all_bids = []
            self.all_vids = []
            self.isb_on_convhull = []
            self.isv_on_convhull = []
            self.v_norms = []
            self.adj_list = []

    def boywerWatson(self):
        d = 100
        far_points = np.array([[0, 0, d], 
                               [d, 0, -d], 
                               [-d, d, -d], 
                               [-d, -d, -d]]) # Only works for R3
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
        
        # Copy the information in db to this instance
        self.delaunay = current_simplices
        self.v2b = db.v2b[:-4]
        all_bids = []
        for s in current_simplices:
            all_bids += s.bids
        # Counting sort
        bid_max = np.max(all_bids)
        self.bid_max = bid_max
        self.b2v = [[] for i in range(bid_max+1)]
        # self.b2s = [[] for i in range(bid_max+1)]
        self.all_bids = []
        self.convhull.isb_on_convhull =[False for i in range(bid_max+1)]
        bid_count = np.zeros(bid_max+1)
        for bid in all_bids:
            bid_count[bid] += 1
        for bid in range(bid_max+1):
            if bid_count[bid] > 0:
                self.b2v[bid] = db.b2v[bid]
                self.all_bids.append(bid)
                # Convex hull boundaries
                if bid_count[bid] == 1:
                    self.convhull.all_bids.append(bid)
                    self.convhull.isb_on_convhull[bid] = True

        # for sid,s in enumerate(current_simplices):
        #     for bid in s.bids:
        #         self.b2s.append(sid)

        return
    
    def processConvHull(self):
        # Convex hull boundary adjacency list
        self.convhull.all_vids = []
        self.convhull.isv_on_convhull = [False for i in range(self.num_v)]
        self.convhull.adj_list = [[] for i in range(self.bid_max+1)]
        self.b2e = [[] for i in range(self.bid_max+1)]
        for bid in self.convhull.all_bids:
            vids = self.b2v[bid]
            for vid in vids:
                self.convhull.isv_on_convhull[vid]=True
            vid_pairs = [[vids[1],vids[0]],
                        [vids[2],vids[1]],
                        [vids[0],vids[2]]]
            for p in vid_pairs:
                eid = self.existEdge(p)
                if eid < 0:
                    eid = self.addEdge(p)
                    self.e2b.append([bid])
                    self.b2e[bid].append(eid)
                else:
                    self.e2b[eid].append(bid)
                    self.b2e[bid].append(eid)

        for bids in self.e2b:
            self.convhull.adj_list[bids[0]].append(bids[1])
            self.convhull.adj_list[bids[1]].append(bids[0])

        for i in range(self.num_v):
            if self.convhull.isv_on_convhull[i]:
                self.convhull.all_vids.append(i)

        self.convhull.v_norms = np.zeros((self.num_v,3))

        for vid in self.convhull.all_vids:
            vids_neighbor = []
            for eid in self.v2e[vid]:
                if self.e2v[eid][0] == vid:
                    vids_neighbor.append(self.e2v[eid][1])
                else:
                    vids_neighbor.append(self.e2v[eid][0])
            vcoords_neighbors = self.vcoords[vids_neighbor,:]-self.vcoords[vid]
            vcoords_neighbors = (vcoords_neighbors.T/np.linalg.norm(vcoords_neighbors.T,axis=0)).T

            for v in vcoords_neighbors:
                self.convhull.v_norms[vid] += v * np.min(vcoords_neighbors @ v.T)
            self.convhull.v_norms[vid] /= np.linalg.norm(self.convhull.v_norms[vid])
            # cov_neighbors = vcoords_neighbors.T @ vcoords_neighbors
            # w,vec = np.linalg.eigh(cov_neighbors)
            # print("vid: {}, neigbor:{}, w: {}".format(vid,vids_neighbor,w))
            # self.convhull.v_norms[vid] = vec[:,np.argmin(w)]

            # avg_neighbors = np.mean(vcoords_neighbors,axis=0)
            # if self.convhull.v_norms[vid]@avg_neighbors>0:
            #     self.convhull.v_norms[vid] *= -1

            # import matplotlib.pyplot as plt
            # from mpl_toolkits.mplot3d.art3d import Poly3DCollection
            # ax = plt.figure().add_subplot(projection='3d')
            # ax.scatter3D(vcoords_neighbors[:,0],vcoords_neighbors[:,1],vcoords_neighbors[:,2])
            # ax.scatter3D(0,0,0,color = 'red')
            # for v in vcoords_neighbors: ax.plot3D([0,v[0]],[0,v[1]],[0,v[2]],color = 'tab:blue')
            # ax.plot3D([0,self.convhull.v_norms[vid,0]],[0,self.convhull.v_norms[vid,1]],[0,self.convhull.v_norms[vid,2]],color='red')
            # ax.plot3D([0,avg_neighbors[0]],[0,avg_neighbors[1]],[0,avg_neighbors[2]],color='green')
            # ax.axis('equal')   
            # plt.show()
            # 1



