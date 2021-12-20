/** @page createMesh
\par createMesh
This is a simple mesh that will load an input file, read the mesh database, and create the desired
mesh.
See \ref AMP::Mesh::Mesh "Mesh" and \ref AMP::Mesh::BoxMesh "BoxMesh" for more info on the mesh
interfaces and options.
To run the test compile the example and then run the command "createMesh inputFile".
If AMP was compiled with silo, the mesh will be written to a silo file named
"output/inputFile.silo".

\par createMesh.cc
\include createMesh.cc

\par Example 1:
The example input shows the creation of a simple cube mesh in 3d.
\include input_createMesh-cube

\par Example 2:
The example input shows the creation of several different mesh types in 2d
\include input_createMesh-2d

\par Example 3:
The example input shows the creation of several different mesh types in 3d
\include input_createMesh-3d

\par Example 4:
The example input shows the creation of a fuel assembly with 8 pins
\include input_createMesh-Multimesh

*/
