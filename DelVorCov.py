from Geometry import *


class DelVorCov(GeometryDatabase):
    def __init__(self, vcoords):
        super().__init__(vcoords)
        self.delaunay = [] # list of simplices
        # self.b2s = [] # boundaries to simplices
        self.convhull = self.ConvexHull() 
        self.vordiag = self.Voronoi() 
        self.boywerWatson()
        self.processConvHull()
        self.processVoronoi()

    class ConvexHull():
        def __init__(self):
            self.all_bids = []
            self.all_vids = []
            self.isb_on_convhull = []
            self.isv_on_convhull = []
            self.v_norms = []
            self.adj_list = []

    class Voronoi():
        def __init__(self) -> None:
            self.vor_coords = []
            self.v2vor = []

    def boywerWatson(self):
        far_points = np.array([[1,0,-np.sqrt(2)], 
                               [-1,0,-np.sqrt(2)], 
                               [0,1,np.sqrt(2)], 
                               [0,-1,np.sqrt(2)]]) * 3 * np.max(np.linalg.norm(self.vcoords,axis = 1))
        print(far_points)
        vcoords = self.vcoords
        far_ids = [len(vcoords),len(vcoords)+1,len(vcoords)+2,len(vcoords)+3]
        vcoords = np.block([[vcoords],[far_points]])
        db = GeometryDatabase(vcoords)
        current_simplices = [Simplex(far_ids, db)]
        for vid in range(len(vcoords)-4):
            # print("vid:{}".format(vid))
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
                    # print("    bid:{}".format(bid))
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

            # Method1
            n1 = np.zeros(3)
            vid1 = vids_neighbor[0]
            vids_neighbor_sort = [vid1]
            vids_neighbor_ = vids_neighbor[1:]

            while len(vids_neighbor_)>0:
                for i,vid2 in enumerate(vids_neighbor_):
                    if self.existEdge([vid1,vid2]):
                        vids_neighbor_sort.append(vid2)
                        vids_neighbor_.pop(i)
                        vid1 = vid2
                        break
            vids_neighbor_sort.append(vids_neighbor[0])
            for i,vid1 in enumerate(vids_neighbor_sort[:-1]):
                vid2 = vids_neighbor_sort[i+1]
                v1 = self.vcoords[vid1,:]-self.vcoords[vid]
                v2 = self.vcoords[vid2,:]-self.vcoords[vid]
                v1 /= np.linalg.norm(v1)
                v2 /= np.linalg.norm(v2)
                n1 -= v1 * np.linalg.norm(np.cross(v1,v2))
            n1 /= np.linalg.norm(n1)

        
            # Method 2

            vcoords_neighbors = self.vcoords[vids_neighbor,:]-self.vcoords[vid]
            vcoords_neighbors = (vcoords_neighbors.T/np.linalg.norm(vcoords_neighbors.T,axis=0)).T

            cov_neighbors = vcoords_neighbors.T @ vcoords_neighbors
            w,vec = np.linalg.eigh(cov_neighbors)
            n2 = vec[:,np.argmin(w)]

            avg_neighbors = np.mean(vcoords_neighbors,axis=0)
            if n2@avg_neighbors>0:
                n2 *= -1

            if np.min(np.arccos(vcoords_neighbors @ n1.T)) > np.min(np.arccos(vcoords_neighbors @ n2.T)):
                self.convhull.v_norms[vid] = n1
            else:
                self.convhull.v_norms[vid] = n2
    
    def processVoronoi(self):
        self.vordiag.vor_coords = np.zeros((len(self.delaunay),3))
        self.vordiag.v2vor = [[] for i in range(self.num_v)]
        for sid,s in enumerate(self.delaunay):
            self.vordiag.vor_coords[sid] = s.cc
            for vid in s.vids:
                self.vordiag.v2vor[vid].append(sid)

    def Amenta(self):
        F_ids = []
        for vid in range(self.num_v):
            if self.convhull.isv_on_convhull[vid]:
                n_plus = self.convhull.v_norms[vid]
            else:
                d = 0
                p_plus = np.zeros(3)
                p_plus_id = -1
                for vorid in self.vordiag.v2vor[vid]:
                    pp = self.vordiag.vor_coords[vorid] 
                    dd = np.linalg.norm(pp - self.vcoords[vid])
                    if dd>d:
                        d = dd
                        p_plus_id = vorid
                F_ids.append(p_plus_id)
                n_plus = self.vordiag.vor_coords[p_plus_id] - self.vcoords[vid]

            d = 0
            p_minus_id = -1
            for vorid in self.vordiag.v2vor[vid]:
                pp = self.vordiag.vor_coords[vorid]
                dd = (pp - self.vcoords[vid]) @ n_plus
                if dd < d:
                    d = dd
                    p_minus_id = vorid
            F_ids.append(p_minus_id)

        
        Fid_max = np.max(F_ids)
        Fcount = [False for i in range(Fid_max+1)]
        for fid in F_ids:
            Fcount[fid] = True
        F_ids = []
        for fid in range(Fid_max+1):
            if Fcount[fid]:
                F_ids.append(fid)
        F = self.vordiag.vor_coords[F_ids]
        voords_new = np.block([[self.vcoords],[F]])
        print(voords_new.shape)
        db = DelVorCov(voords_new)
        crust = []
        for bid in db.all_bids:
            is_crust = True
            for vid in db.b2v[bid]:
                if vid > self.num_v-1:
                    is_crust = False
                    break
            if is_crust:
                crust.append(bid)
        return db, crust


                    
