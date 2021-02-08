## Escaping Villa Pisani

Matlab codes for the simulation of the Patlak-Keller-Segel equations inside Villa Pisani and for circular target, with learning dynamics (mean field optimal escape/search with energy, time and collision costs).

Related publication:

Pezzotta, A., Adorisio, M., & Celani, A. "Chemotaxis emerges as the optimal solution to cooperative search games". *Phys. Rev. E*, **98**, 42401 (2018)
<!-- https://doi.org/10.1103/PhysRevE.98.042401 -->


### Run

Script uses PDETools package

```bash
matlab PDE_maze
```

Working on Matlab 2018b

Comment/uncomment relevant lines in order to change the geometry.

### Other files

Files `centered_*.dat` contain the 2D coordinates of the boundaries of the *Villa Pisani* domain. The domain is not simply connected, so `centered_1.dat` contains the coordinates of the outer domain. The other files contain the coodinates of the domains to *subtract* to the outer one.
