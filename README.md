# Electrostatic_phase_plate
Source code for the image analysis of electron holography characterization of an electrostatic phase plate

**In the repository**:
Two Python Files for Image analysis:
 * "_Algo1" is the first stitching, where the raw data is first analyzed by means of phase reconstruction and phase unwrapping
 * "_Algo2" is the second image analysis, where the slope in phase is first measured and an average slope is then estimated. The average slope is then used to realign the holograms, where the average phase of each hologram is first subtracted and then the overall trend is added to that hologram. This is repeated for the entirity of the images produced from "_Algo1".

**The packages/libraries used for analysis are:**
 * Numpy (NP)
 * Matplotlib for plotting and image show
 * Hyperspy (Holography, https://hyperspy.readthedocs.io/en/stable/user_guide/electron_holography.html)
 * Skimage for phase unwrapping (Hyperspy's Phase Unwrapping is that from Skimage)

Contact : madsslarsen95@gmail.com
