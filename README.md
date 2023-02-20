# Poisson-Equation-solved-with-DEAL.II
This Repository contains code for solving Poisson's Equations in three geometries of Square, Circle and a Ring. 


## Poisson's Equation 
‍‍
```math
-\Delta u = f \;\; in \;\; \Omega
```

```math
u = 0 \;\; in \;\; \partial \Omega
```

## Compiling
for compling the *.cc codes following commands must be executed:

```
sudo cmake -DDEAL_II_DIR="path\to\DEAL.II_directory" .
```
then:
```
make
```
the executable program will be compipled and after running it the **solution.vtk**, **solution.gpl** and **grid.svg** will be generated.

## GnuPlot results
for visualizing the 3D surface the follwing commands must be executed in GnuPlot shell:
```
gnuplot
set style data line
splot "solution.gpl"
```

## Results
### Cube
![cube grid](https://github.com/DaraSamii/Poisson-Equation-solved-with-DEAL.II/blob/master/images/grid_cube.png)
![cube contour](https://github.com/DaraSamii/Poisson-Equation-solved-with-DEAL.II/blob/master/images/cube_contour.png)
### Circle
![circle grid](https://github.com/DaraSamii/Poisson-Equation-solved-with-DEAL.II/blob/master/images/circle_grid.png)
![circle contour](https://github.com/DaraSamii/Poisson-Equation-solved-with-DEAL.II/blob/master/images/circle_contour.png)
### Ring
![ring grid](https://github.com/DaraSamii/Poisson-Equation-solved-with-DEAL.II/blob/master/images/ring_grind.png)
![ring contour](https://github.com/DaraSamii/Poisson-Equation-solved-with-DEAL.II/blob/master/images/ring_contour.png)
