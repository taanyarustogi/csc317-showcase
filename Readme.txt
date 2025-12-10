Title: The Jelly Cubes

Personal Information:
Full Name: Tanya Rustogi
UtorID: rustogit
Student Number: 1007828665
Assignment Augmented: A8

Instructions
1. Compilation:
    cd code
    mkdir build
    cd build
    cmake ..

2.  Build: 
    Open the .sln file in the build folder and press build to make sure it's built.

3.  Run: 
    navigate to the build folder or the build/debug folder
    ./masssprings_dense ../jelly2.json

4.  Interaction:
    Press 'A' to start the animation loop.
    Press 'O' repeatedly to sequentially release the cubes

Description:
This Piece edits the Mass-Spring System Assignment to simulate three jelly cubes that interact together to create a jelly cube tower. 

Features Added:

1. Key-Press Release
    - The cubes are initially fixed in place. Pressing the 'O' key releases the next available cube in the stack.
    - Technical implementation: a state tracking vector (std::vector<bool> cube_is_fixed) and conditional application of the external gravity force.
    - Code Location: main.cpp

2. Floor Collision
    - To ensure that the cubes collide with a floor upon release.
    - Technical implementation: any vertex below 0 is adjusted to be above the floor and the reverse velocity simulates a bounce.
    - Code Location: fast_mass_springs_step_dense.cpp

4. Vertical Cube-to-Cube Collision Resolution
    - Collision logic that ensures that the cubes stack on top of each other.
    - Technical implementation: any vertex that penetrates to the next cube is increased to stop overlap and create a bounce.
    - Code Location: fast_mass_springs_step_dense.cpp

5. Custom Geometry
    - Created a custom geometry file (jelly-2.json and jelly2.obj), defining the springs and vertices for 3 cubes.
