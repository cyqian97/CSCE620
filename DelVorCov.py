from Geometry import *


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
        
        # Copy the information in db to this instance
        self.v2b = db.v2b[:-4]
        all_bids = []
        for s in current_simplices:
            all_bids += s.bids
        # Counting sort
        bid_max = np.max(all_bids)
        self.b2v = [[] for i in range(bid_max+1)]
        self.all_bids = []
        bid_count = np.zeros(bid_max+1)
        for bid in all_bids:
            bid_count[bid] += 1
        for bid in range(bid_max+1):
            if bid_count[bid] > 0:
                self.b2v[bid] = db.b2v[bid]
                self.all_bids.append(bid)

        self.delaunay = current_simplices

        return