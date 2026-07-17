import numpy as np

class Polymer():

    def __init__(self, nMonomers, bondLength, m, k, box, kBT, seed):
        self.nMonomers = nMonomers
        self.bondLength = bondLength
        self.m = m
        self.k = k
        

        self.box = box
        self.kBT = kBT

        self.rng = np.random.default_rng(seed)

        self.generateLinearPositions()
        self.v = np.zeros((nMonomers, 3))
        self.f = np.zeros((nMonomers, 3))

    def generatePolymerPositions(self):
        res = []
        self.r[0] = self.rng.uniform(0, self.box, size=(3))
        for i in range(1, self.nMonomers-1):
            self.r

    def generateLinearPositions(self):

        self.r = np.zeros((self.nMonomers, 3))
        self.r[0] = self.rng.uniform(0, self.box, size=3)

        theta, alpha = self.rng.uniform(0, 2*np.pi, size=2)
        direction = np.array([np.sin(theta)*np.cos(alpha), np.sin(theta)*np.sin(alpha), np.cos(theta)])
        step = self.bondLength * direction

        for i in range(1, self.nMonomers):
            self.r[i] = self.r[i-1]+step*self.rng.uniform(0.95,1.05)
            self.r %= self.box

    def vecDiffInBC(self, idx1, idx2):
        vecNoBC = self.r[idx2] - self.r[idx1]
        return (vecNoBC - self.box*np.round(vecNoBC/self.box))
    
    def distInBC(self, idx1, idx2):
        vec = self.vecDiffInBC(idx1, idx2)
        return np.sqrt(np.sum(vec**2))