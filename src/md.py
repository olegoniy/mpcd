import numpy as np

def verlet_step(polymer, dt, force_func):
    polymer.v += dt*polymer.f/polymer.m

    polymer.r += dt * polymer.v
    polymer.r %= polymer.box

    force_func(polymer)

