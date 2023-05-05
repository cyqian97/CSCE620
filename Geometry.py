import numpy as np
import itertools

class GeometryDatabase:
    def __init__(self, vcoords):
        self.vcoords = vcoords
        self.num_v = vcoords.__len__()
        self.bid_max = 0
        self.eid_max = 0
        self.v2b = [[] for i in range (self.num_v)] # List of list of boundaries attached to each vertex
        self.v2e = [[] for i in range (self.num_v)]
        self.b2v = [] # List of list of vetices of each boundary
        self.b2e = []
        self.e2v = []
        self.e2b = []
        self.all_bids = []
        self.all_eids = []

    def getVCoords(self,vids):
        return np.array([self.vcoords[v] for v in vids])
    
    def getBoundaryVCoords(self,bid):
        return self.getVCoords(self.b2v[bid])
    
    def existBoundary(self,vids):
        b = set.intersection(*[set(self.v2b[id]) for id in vids])
        if len(b) == 0:
            return -1
        elif len(b) == 1:
            return list(b)[0]
        else:
            raise Exception("Multiple boundaries are found: {}".format(b))

    def addBoundary(self,vids):
        if self.existBoundary(vids) < 0:
            self.b2v.append(vids)
            bid = len(self.b2v)-1
            self.bid_max = bid
            self.all_bids.append(bid)
            for v in vids:
                self.v2b[v] = self.v2b[v] + [bid]
            return bid
        else:
            return -1
    
    def existEdge(self,vids):
        if not len(vids)==2:
            raise Exception("Need 2 points but given {}".format(len(vids)))
        b = set.intersection(set(self.v2e[vids[0]]),set(self.v2e[vids[1]]))
        if len(b) == 0:
            return -1
        elif len(b) == 1:
            return list(b)[0]
        else:
            raise Exception("Multiple boundaries are found: {}".format(b))

    def addEdge(self,vids):
        if self.existEdge(vids) < 0:
            self.e2v.append(vids)
            eid = len(self.e2v)-1
            self.eid_max = eid
            self.all_eids.append(eid)
            for v in vids:
                self.v2e[v].append(eid)
            return eid
        else:
            return -1
        
    def faceNormal(self,bid):
        vcoords = self.getBoundaryVCoords(bid)
        n = np.cross(vcoords[2]-vcoords[1],vcoords[1]-vcoords[0])
        return n / np.linalg.norm(n)


class Simplex:
    def __init__(self, vids, db):
        """Initialize Simplex insstance

        Args:
            points (np.array(dim+1,dim)): vertices of the simplex
        """
        vertices = db.getVCoords(vids)

        if vertices.shape[0]-vertices.shape[1] == 1:
            self.dim = vertices.__len__()-1
            self.vids = vids
            self.vertices = vertices
            self.bids = []
        else:
            raise Exception("The dimension of the vertices is incorrect. Current shape: ({},{})".format(vertices.shape(0),vertices.shape(1)))
        
        for c in itertools.combinations(vids,self.dim):
            bid = db.existBoundary(c)
            if bid < 0:
                self.bids.append(db.addBoundary(c))
            else:
                self.bids.append(bid)
        self.cc, self.cr = self.calcCircumcircle() # Circumcenter and circumsphere radius



    def cayleyMengerMatrix(self):
        m = np.zeros((self.dim+2,self.dim+2))
        m[1:,0] = 1
        m[0,1:] = 1
        for i in range(self.dim+1):
            for j in range(self.dim+1):
                m[i+1,j+1] = np.linalg.norm(self.vertices[i]-self.vertices[j])**2
        return m


    def calcCircumcircle(self):
        """Calculate the center and the radius of the circumcircle

        Returns:
            cc (np.array): _description_
        """
        m = self.cayleyMengerMatrix()
        q = np.linalg.inv(m)
        qq = q[0,1:]
        cc = qq @ self.vertices
        cr = np.sqrt(-q[0,0]/2)
        return cc, cr
    
    def isInCircumcicle(self,p):
        return np.linalg.norm(p-self.cc)<self.cr+1e-30
    
# s = Simplex(np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1]]))
# c,r = s.calcCircumcircle()
# print(c)
# print(r)