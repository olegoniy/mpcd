import sys
sys.path.append("../src")

from system import System
from dynamics import streaming, distributeToCells, collision

system = System(
    N = 50, 
    box = [3, 3, 3], 
    a = 1.0, 
    h = 0.1, 
    m = 1.0,  
    kBT = 1.0, 
    alpha_deg=30, 
    seed=12345
    )


#print(system.r)
collision(system)
