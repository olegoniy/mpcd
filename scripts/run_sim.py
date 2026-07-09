import sys
sys.path.append("../src")
print(sys.path)

from system import System
from dynamics import streaming
from observables import system_impuls, system_kinetic

system = System(
    N = 10000, 
    box = [5, 5, 5], 
    a = 1.0, 
    h = 0.1, 
    m = 1.0, 
    kBT = 1.0, 
    alpha_deg=30, 
    seed=12345
    )

for j in range(100):
    streaming(system)
    print(f"Iterration {j+1}|\nPosition: {system.r[0]}\nImpuls: {system_impuls(system)} \nKinetic energy: {system_kinetic(system)}\n"
    )


