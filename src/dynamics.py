import numpy as np

def streaming(system):
    system.r += system.h * system.v
    system.r %= system.box
