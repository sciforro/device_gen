# device_gen
This is a code repository collecting utilities and functions for creating devices for microfabrication. It focuses on mostly micro-electroed-arrays for which automated routing is under continued development.

## Tutorial #1: Basics of Shapely and Phidl:
Two Python packages that synergize very well for creating device layouts are 
- Shapely: low-level, extremely powerful geometry library
- Phidl : high-level, CAD-oriented package with a focus on photonics, with a large number of quality-of-life functionalities and a gds-ii format exporter)

Video tutorial about the basics is linked below

[![Phidl and Shapely basics](https://img.youtube.com/vi/8cDB7dCHEBI/0.jpg)](https://www.youtube.com/watch?v=8cDB7dCHEBI) 

with the juyter notebook to follow along.

## Tutorial #2: Deep dive into auto-routing:
The following tutorials goes into the philosophy and algorithms behind automated routing. Given an input device, the shape and location of both electrode and pads, the RoutingLayout class will try to automatically route pads and electrodes via creating an implicit graph in the device and use maximal flow algorithms on the graph to compute node-disjoint paths (non-intersecting paths between electrodes and pads, which is a key requirement when routing electronics as metal interconnects are on the same metal layer and cannot intersect). If you are interested in how to use the routing framework, skip this video as here I derive the principles and functions to achieve automated routing.

[![Automated routing](https://img.youtube.com/vi/23VDh5naxl4/0.jpg)](https://www.youtube.com/watch?v=23VDh5naxl4) 

You can follow along in this jupyter notebook

## Tutorial #3: Device example 1:
Here I reproduce a device by McDonald et. al published at http://dx.doi.org/10.1016/j.bios.2023.115223 
The ease of design and automated routing is showcased. 

[![Example device - 1](https://img.youtube.com/vi/I84Mfm7iDBs/0.jpg)](https://www.youtube.com/watch?v=I84Mfm7iDBs)



