## Project: Kinematics Pick & Place

### Project overview

#### What is this?
This is my implementation of the Udacity Robotics Nanodgree Pick-and-Place project.

#### Problem statement
A six degree of freedom (six jointed) robot arm (Kuka kr210) is in front of a shelf with 9 compartments. Cylinders spawn one at a time on one compartment at random. What needs to be done is for us to determine the correct angles and positions for each of the robot arm's joints so that the joints (and the robot as a whole) can trace a path to grab the spawned cylinder and place it into a nearby bucket. A gif of what needs to happen is shown below(courtesy of Udacity Robotics Nanodegree material):

![what_needs_to_happen](./writeup_images/what_needs_to_happen.gif)

#### Solution and files of note
- The writeup of this project which includes the problem solving steps can be found here: [WRITEUP](./pick_and_place_writeup.md).

- A video of the robot arm in action can be found here: [PICK AND PLACE VIDEO](./video/kuka_pick_and_place.ogg)

- The solution code can be found here: [IK_server.py](./IK_server.py). It is necessary to have Gazebo, ROS, and MoveIt! to be all available to run this code. Running this code will produce an `error_list.json` file that can be used later on to plot error rates.
    - How to run the code:
        1. Install ROS Kinetic Kame on your Ubuntu machine
        2. Follow the instructions in the [original repo's readme](https://github.com/udacity/RoboND-Kinematics-Project).
        3. Look for the `inverse_kinematics.launch` file under `/RoboND-Kinematics-Project/kuka_arm/launch`. Make sure to set the `demo` flag to `false`.
        4. Overwrite `IK_server.py` in `~/catkin_ws/src/RoboND-Kinematics-Project/kuka_arm/scripts` with the above mentioned [IK_server.py](./IK_server.py) file.
        5. Open a terminal and launch the project by invoking the following commands
            ```bash
            $ cd ~/catkin_ws/src/RoboND-Kinematics-Project/kuka_arm/scripts
            $ ./safe_spawner.sh 
            ```
        6. Open a separate terminal and run the script we copy and pasted earlier:
            ```bash
            $ cd ~/catkin_ws/src/RoboND-Kinematics-Project/kuka_arm/scripts
            $ rosrun kuka_arm IK_server.py   
            ```
        7. This should produce an `error_list.json` file, which lists the error rates.  If desired, an error histogram can be produced by running the `error_plotter.py` script which is described below.

- Stand alone solution code with sample data
    - This script is pretty much the same as the aforementioned solution code.  However, this can be run without having to open ROS, Gazebo, or MoveIt! as long as its dependencies which include numpy, sympy, and the sample data are available. This does not produce an `error_list.json` file, however it prints out error rates to the console.
    - [IK.py - stand alone script](./IK.py)
    - [sample_data.py - necessary to run IK.py](./sample_data.py)

- Error plotter script
    - This script will generate the error plot histograms that were used in this writeup. It requires numpy as well as [plotly](https://plot.ly/) to run. It will also need an `error_list.json` file that can be obtained by running the aforementioned [IK_server.py](./IK_server.py) solution with ROS, Gazebo, and MoveIt!. To run the error plotter script that is found below, simply call it with python. Make sure that error_list.json is in the same directory.
    ```python
    $ python ./error_plotter.py
    ```
    - [error_plotter.py](./error_plotter.py)
    - [error_list.json](./error_list.json)