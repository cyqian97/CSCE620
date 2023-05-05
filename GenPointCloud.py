import numpy as np
from time import time

all_types = ["cat","rand","randtorus","diamond","prism","icosahedron"]

def genPointCloud(type,n=0):
    if type == "cat":
        vcoords = np.loadtxt("cat.txt")
        vcoords -= np.mean(vcoords)
        vcoords /= np.std(vcoords)
        selected = np.random.randint(0,vcoords.shape[0],(n-1,))
        
        return np.array([vcoords[i] for i in set(selected)])
    elif type == "rand":
        seed = int(time())
        np.random.seed(seed)
        print("Random seed: {}".format(seed))
        # np.random.seed(1683256853)
        vcoords = np.random.uniform(0,1,(n,3))
    elif type == "torus":
        
        seed = int(time())
        np.random.seed(seed)
        print("Random seed: {}".format(seed))
        # np.random.seed(1683252544)
        r1 = 0.5
        r2 = 0.25
        n1 = int(np.sqrt(n)/np.sqrt(2))
        n2 = int(n/n1)
        angles1 = np.linspace(0,np.pi*2*(n1-1)/n1,n1)
        angles2 = np.linspace(0,np.pi*2*(n2-1)/n2,n2)
        noise = 0.02
        vcoords = []
        i = 0
        for a in angles2:
            for b in angles1:
                vcoords.append([
                    ((r1+np.random.uniform(-noise, noise))+ (r2+np.random.uniform(-noise, noise)) * np.cos(b+np.random.uniform(-noise, noise))) * np.cos(a+np.random.uniform(-noise, noise)),
                    ((r1+np.random.uniform(-noise, noise))+ (r2+np.random.uniform(-noise, noise)) * np.cos(b+np.random.uniform(-noise, noise))) * np.sin(a+np.random.uniform(-noise, noise)),
                    (r2+np.random.uniform(-noise, noise)) * np.sin(b+np.random.uniform(-noise, noise))
                ])
                i += 1
        return np.array(vcoords)

    elif type == "randtorus":
        r1 = 0.5
        r2 = 0.25
        seed = int(time())
        np.random.seed(seed)
        print("Random seed: {}".format(seed))
        angles = np.random.uniform(0,2*np.pi,(n,2))
        vcoords = np.block([
        [(r1+ r2 * np.cos(angles[:,1])) * np.cos(angles[:,0])],
        [(r1+ r2 * np.cos(angles[:,1])) * np.sin(angles[:,0])],
        [r2 * np.sin(angles[:,1])]
        ]).T
    elif type == "diamond":
        vcoords = np.array([[0,0,1],
                            [np.cos(np.pi/6),np.sin(np.pi/6),0],
                            [np.cos(5*np.pi/6),np.sin(5*np.pi/6),0],
                            [0,-1,0],
                            [0,0,-1]])*1.5
    elif type == "prism":
        vcoords = np.array([[np.cos(np.pi/6), np.sin(np.pi/6), 1],
                            [np.cos(5*np.pi/6), np.sin(5*np.pi/6), 1],
                            [0, -1, 1],
                            [np.cos(np.pi/6), np.sin(np.pi/6), -1],
                            [np.cos(5*np.pi/6), np.sin(5*np.pi/6), -1],
                            [0, -1, -1]])
    elif type == "icosahedron":
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
    return vcoords
