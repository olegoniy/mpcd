import numpy as np


class System:
    def __init__(self, N, box, a=1.0, h=0.1, m=1.0, kBT=1.0, alpha_deg=228, seed=69):
        self.N = N
        self.box = np.array(box, dtype=float)
        self.a = a
        self.h = h
        self.m = m
        self.kBT = kBT
        self.alpha = np.deg2rad(alpha_deg)

        self.rng = np.random.default_rng(seed) # to generate random numbers

        self.r = self.rng.uniform(0, self.box, size=(N, 3))
        self.v = self.rng.normal(0, np.sqrt(kBT / m), size=(N, 3))

        self.v -= self.v.mean(axis=0) #to get drift away it seems, it will help with temperature acurracy

        