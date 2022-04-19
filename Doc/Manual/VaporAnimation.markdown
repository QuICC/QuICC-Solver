Creating Animations in VAPOR {#pManVaporAnim}
============================

This document is a reference guide for creating animations using visualizations in VAPOR (*Visualization and Analysis Platform for Ocean, Atmosphere, and Solar Researchers*) open source software from the National Center for Atmospheric Research's Computational and Information Systems Lab. Further, the visualizations are captured in a movie using *FFMPEG*, A complete, cross-platform solution to record, convert and stream audio and video.

Introduction
------------

An animation with VAPOR works much like a cartoon, meaning it will require multiple visualization state files statistically close together (similar). This means that the system being visualized must have reached a steady state solution, so a short simulation producing a large number of state files will be sufficient for the animation. Using these state files, visualization state files are produced using a shell, then the files are visualized in VAPOR as normal, but with slight modifications. Using VAPOR's capture tool, a sequence of images can be captured and finally put into an animation using ffmpeg. 

Note, to make an animation, the system must have reached a steady state, else we will be dealing with states that are so statistically unstable that each sequential image will be drastically different. This will cause problems with VAPOR's transfer function and capture tools. Besides that, the information before reaching a steady state is not particularly illuminating. 

Generating Visualization Files
------------------------------

The next step is to generate visualization files (vaporVisState.hdf5) from the state files (state.hdf5). This process is the same as presented in refman.pdf, but must be done individually for each state file. To save time, a shell script has been implemented to step through this process: 


*SEE SLURM FILES:*

- visStates.slurm
- vaporVisStates.slurm

When this is finished, there will be *n* visualization state files (vaporVisState000n.hdf5), and three data files including information about the *x,y* and *z* axis for visualization. At this point, all files must be copied to a local directory to prepare for visualization. 

Importing Files Together in VAPOR
---------------------------------

Once all visualization files (vaporVisState000n.hdf5, visState0000_vapor_(x,y,z).dat) are copied into a local directory, the same instructions presented in refman.pdf can be followed with a few slight modifications. 


Create the VDF file:

- **TFF scheme:**
#$ > ncdfvdfcreate -periodic 1:1:0 -level 3 -xcoords visState0000_vapor_z.dat -ycoords vis-
State0000_vapor_y.dat -zcoords visState0000_vapor_x.dat vaporVisState0000.hdf5 vaporVisState0001.hdf5 \dots vaporVisState00n.hdf5 vis_me.-
vdf

- **Populate NC data:**
#$ > ncdf2vdf -level 3 vaporVisState0000.hdf5 vaporVisState0001.hdf5 \dots vaporVisState00n.hdf5 vis_me.vdf

- **Visualize:**
#$ > vaporgui vis_me.vdf


Now, VAPOR will open with all *n* timesteps in sequential order. Simply use the arrows at the top of the page to step through each image. 

Capturing Images in VAPOR
-------------------------

Now we wish to capture an image for each timestep that will later be `bound together' to make the animation. To capture the images, VAPOR has a nice tool that will take a screenshot everytime the image in the region changes. It will name the images (.jpg or .tiff) in sequential order and save them to a directory of your choice. To do this:

*Capture --> Begin image capture sequence in visualizer No.*

Now choose your directory, name the image, and play the animation. When finished:

*Capture --> End image capture sequence in visualizer No.*

Finally, check that the directory is full of your desired images. 

Creating a Movie with FFMPEG
----------------------------

The last step is to congregate the image files (.jpg or .tiff) into a movie file (.mp4). This documentation uses ffmpeg for this, but note it can be done with different software. The command to make a movie out of multiple images is as follows:

`#$ >  ffmpeg -framerate 8 -i $~$/path2images/image_name\%4d.jpg -c:v libx264 -vf `
` "scale=854:trunc(ow/a/2)*2" -r 30 -pix_fmt yuv420p movie_name.mp4 `

This will save the .mp4 file to the directory in which ffmpeg was executed. For more documentation on ffmpeg see sources below.

Sources
-------

VAPOR Visualization & Analysis Platform
[vapor.ucar.edu](https://www.vapor.ucar.edu/)

FFMPEG
[ffmpeg.org](https://www.ffmpeg.org/documentation.html)

\author Derek Driggs
\author Louie Long
\date June 29, 2015

