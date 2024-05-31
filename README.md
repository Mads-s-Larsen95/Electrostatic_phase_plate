# Electrostatic_phase_plate
Source code for the image analysis of electron holography characterization of an electrostatic phase plate

**In the repository**:
Python file for image analysis of the electrostatic phase plate in action for 2V applied. The inputs are .dm4 (Gatan) objects. 
 * Creates the 2D map observed in the main article
   *   The location in (x,y) is found from the source of the image
   *   To measure the offset, the last 30 rows/columns is used as an offset for the neighboring image
   *   Whenever the y position is changed (rows), the offset is taken from the last 30 rows of the first image in the adjacent image
   *   Whenever the x position is changed (columns), the offset is taken from the last 30 columns of image coming before


**The packages/libraries used for analysis are:**
 * Numpy (NP)
 * Matplotlib for plotting and image show
 * Hyperspy (Holography, https://hyperspy.readthedocs.io/en/stable/user_guide/electron_holography.html)
 * Skimage for phase unwrapping (Hyperspy's Phase Unwrapping is that from Skimage)

Contact : madsslarsen95@gmail.com
