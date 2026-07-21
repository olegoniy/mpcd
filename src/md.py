import numpy as np

def verlet_step(polymer, dt, force_func):
    polymer.v += dt*polymer.f/polymer.m * 0.5

    polymer.r += dt * polymer.v
    polymer.r %= polymer.box

    force_func(polymer)
    polymer.v += dt*polymer.f/polymer.m * 0.5

