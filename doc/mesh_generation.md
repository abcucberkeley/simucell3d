# Mesh generation

Generating meshes is a cumbersome task. However, there are some tools that can facilitate this process. 

## Paraview
The first one is the [Paraview](https://www.paraview.org/download/) software.


You can select a cell from a tissue by first clicking on the `Ctrl` key and then clicking on the `s` key and finally by making a left click on the cell you want to extract:
<p align="center">
    <img src="./img/cell_selection.png">
</p>

You can then extract this cell from the rest of the mesh by using the filter `extract selection`:
<p align="center">
    <img src="./img/extract _from_selection.png">
</p>
The extracted cell can then be translated, rotated and scaled to the desired position with the `transform` filter:
<p align="center">
    <img src="./img/transform_filter.png">
</p>

Different meshes can be regrouped in one mesh with the `append datasets` filter:
<p align="center">
    <img src="./img/append_filter.png">
</p>

<br>

## MeshLab
The meshlab software can also be used to generate meshes. It is available [here](http://www.meshlab.net/).
<p align="center">
    <img src="./img/meshlab.png">
</p>

It contains various filters that can be used to smooth the meshes, remove points or triangles or regenerate the mesh surfaces.
<p align="center">
    <img src="./img/meshlab_filters.png">
</p>
<br>

## Blender

If you have to generate a mesh from start Blender is probably the best option. It is a free and open source 3D creation suite designed for videos and games production. It is available [here](https://www.blender.org/download/).

Blender offers a wide range of tools to create meshes, to sculpt the meshes and to position their nodes:
<p align="center">
    <img src="./img/blender.png">
</p>
