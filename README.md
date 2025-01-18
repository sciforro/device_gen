# device_gen
This is a code repository collecting utilities and functions for creating devices for microfabrication. It focuses on mostly micro-electroed-arrays for which automated routing is under continued development.

## Tutorial #1: Basics of Shapely and Phidl:
Two Python packages that synergize very well for creating device layouts are 
- Shapely: low-level, extremely powerful geometry library
- Phidl : high-level, CAD-oriented package with a focus on photonics, with a large number of quality-of-life functionalities and a gds-ii format exporter)

Video tutorial about the basics is linked below

[![Phidl and Shapely basics](https://img.youtube.com/vi/8cDB7dCHEBI/0.jpg)](https://www.youtube.com/watch?v=8cDB7dCHEBI) 

with the [jupyter notebook](https://github.com/sciforro/device_gen/blob/main/notebooks_and_code/1_Video_tutorial_basics.ipynb)
 to follow along.

## Tutorial #2: Deep dive into auto-routing:
Make sure to put utilities.py in the same folder as your jupyter notebook

The following tutorials goes into the philosophy and algorithms behind automated routing. Given an input device, the shape and location of both electrode and pads, the RoutingLayout class will try to automatically route pads and electrodes via creating an implicit graph in the device and use maximal flow algorithms on the graph to compute node-disjoint paths (non-intersecting paths between electrodes and pads, which is a key requirement when routing electronics as metal interconnects are on the same metal layer and cannot intersect). If you are interested in how to use the routing framework, skip this video as here I derive the principles and functions to achieve automated routing.

[![Automated routing](https://img.youtube.com/vi/23VDh5naxl4/0.jpg)](https://www.youtube.com/watch?v=23VDh5naxl4) 

You can follow along in this [jupyter notebook](https://github.com/sciforro/device_gen/blob/main/notebooks_and_code/2_Video_tutorial_on_routing.ipynb)

## Tutorial #3: Device example 1:
Make sure to put utilities.py (in the notebooks_and_code folder) in the same folder as your jupyter notebook
Here I reproduce a device by McDonald et. al. published at http://dx.doi.org/10.1016/j.bios.2023.115223 

The ease of design and automated routing is showcased. 

[![Example device - 1](https://img.youtube.com/vi/I84Mfm7iDBs/0.jpg)](https://www.youtube.com/watch?v=I84Mfm7iDBs)

The corresponding [jupyter notebook](https://github.com/sciforro/device_gen/blob/main/notebooks_and_code/3_Video_tutorial_device_McDonald.ipynb) is found here. The final device will look like this: <img align="center" height="450" src="https://github.com/sciforro/device_gen/blob/main/images/dev1.png">


## Tutorial #4: Device example 2:
Make sure to put utilities.py in the same folder as your jupyter notebook
Here I reproduce a device by Park et. al. published at https://www.science.org/doi/10.1126/sciadv.abf9153 

Specifically, I show how to create the serpentine structure shown in Figure S19 and embed it in the entire device. The routing is performed in a similar fashion to that shown in the publication.

[![Example device - 2](https://img.youtube.com/vi/OHwyicW2AnE/0.jpg)](https://www.youtube.com/watch?v=OHwyicW2AnE)

The corresponding [jupyter notebook](https://github.com/sciforro/device_gen/blob/main/notebooks_and_code/4_Video_Tutorial_device_Park.ipynb)  is found here and the final device will look like this: <img align="center" height="450" src="https://github.com/sciforro/device_gen/blob/main/images/dev2.png">

## Tutorial #5: Spiral device creation, simulation and optimization:
Make sure to put utilities.py and fenics_simulation.py  (in the notebooks_and_code folder) in the same folder as your jupyter notebook

Here, I reproduce in broad strokes the spiral device publish by myself and co-author X. Yang at https://www.nature.com/articles/s41587-023-02081-3 

The relevant part of the device is created. Then, it is converted to a mesh for finite-element simulation with FEniCS, a finite-element solver package written in Python. I show how to prepare the mesh to easily import it into FEniCS, and how to mark the relevant parts of the device during the mesh creation so that they can be used automatically to set boundary condition (e.g no slip/rotation at the device edge). Finally, I show how to simulate the device deformation under its own weight (or with the addition of an organoid), and how to process and visualize the simulation results. I show how to interpret the results from a physical stand point in order to guide the further optimization of the device.

[![Device simulation and optimization](https://img.youtube.com/vi/UfntaHuc6hc/0.jpg)](https://www.youtube.com/watch?v=UfntaHuc6hc)

The tutorial [jupyter notebook](https://github.com/sciforro/device_gen/blob/main/notebooks_and_code/5_Video_tutorial_Kirigami_full.ipynb)  is found here. The strain fields overlaid on the deformation are visualized with Blender, using the strain_mapper.blend file in the notebook_and_code folder, and can be explored interactively: 

<p float="left">
<img align="center" height="180" src="https://github.com/sciforro/device_gen/blob/main/images/defo1.png">
<img align="center" height="180" src="https://github.com/sciforro/device_gen/blob/main/images/defo2.png"></p>







