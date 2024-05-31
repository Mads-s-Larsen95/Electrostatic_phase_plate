# -*- coding: utf-8 -*-
"""
Created on Tue Nov  7 08:44:08 2023

@author: maslar
"""
"""
Libraries to use
"""
"""
Functions to use
"""




from skimage.restoration import unwrap_phase
import os
from matplotlib_scalebar.scalebar import ScaleBar
from matplotlib.widgets import Cursor
import natsort
import hyperspy.api as hs  # Loading + Holo package
import numpy as np  # Arrays

import matplotlib.pyplot as plt  # Plotting
import glob  # Load path
def Get_sb_and_Size(im):
    """
    

    Parameters
    ----------
    im : The initial hologram acquired at a certain biprism voltage and rotation, which should be the same for all
    holograms (object + refence) in that particular for loop.

    Returns
    -------
    sb_position : Position of the sb defined by the user (x,y). 
    sb_size : The size of the mask applied for reconstruction - the WPO approximation is used, i.e.
    the radius of the mask is 1/3 the distance from the sideband to the centerband (direct beam)

    """
    plt.close("all")
    fft = im.fft(True)
    fft_plot = fft.data
    rr, cc = np.shape(fft_plot)
    fig, axs = plt.subplots()
    axs.imshow(np.log(abs(fft_plot)), cmap="gray")
    cursor = Cursor(axs, useblit=True, color='k', linewidth=1)
    zoom_ok = False
    print('Zoom or pan to view, \npress spacebar when ready to click:')
    while not zoom_ok:
        zoom_ok = plt.waitforbuttonpress()
    print('\n Select center of sb')
    points = plt.ginput(1)
    points = np.array(points, dtype=int)
    sb_pos1 = points[0]
    # Create ROI with h=roi_w, w=roi_w (square) around the position chosen by user
    roi_w = 50
    fft_ROI = fft_plot[sb_pos1[1]-roi_w:sb_pos1[1] +
                       roi_w, sb_pos1[0]-roi_w:sb_pos1[0]+roi_w]

    # Find maximum in ROI
    # find max position [i,j]
    max_idx = np.unravel_index(np.log(abs(fft_ROI)).argmax(), fft_ROI.shape)
    max_idx = np.asarray(max_idx[::-1])

    sb_pos1 = np.array([sb_pos1[1]-roi_w + max_idx[1],
                       sb_pos1[0]-roi_w+max_idx[0]])

    x_len = int(abs(sb_pos1[1]-cc/2))
    y_len = int(abs(sb_pos1[0]-rr/2))
    L = np.sqrt(x_len**2+y_len**2)
    sb_size = L/3

    if sb_pos1[0] < rr/2 and sb_pos1[1] < cc/2:  # Upper sb chosen to the left
        Txt = "Upper sb to the left"
        sb_position = np.array([y_len, x_len])
    if sb_pos1[0] < rr/2 and sb_pos1[1] > cc/2:  # Upper sb chosen to the right
        Txt = "Upper sb to the left"
        sb_position = np.array([y_len, int(rr-x_len)])
    if sb_pos1[0] > rr/2 and sb_pos1[1] > cc/2:  # Lower sb chosen to the right
        Txt = "Upper sb to the left"
        sb_position = np.array([int(rr-y_len), int(rr-x_len)])
    if sb_pos1[0] > rr/2 and sb_pos1[1] < cc/2:  # Lower sb chosen to the left
        Txt = "Upper sb to the left"
        sb_position = np.array([int(rr-y_len), x_len])
    print("\n", Txt + ", points chosen : ", sb_pos1)
    axs.plot(sb_pos1[1], sb_pos1[0], 'r.')
    circle = plt.Circle((sb_pos1[1], sb_pos1[0]),
                        sb_size, color="r", lw=1.5, fill=False)
    axs.add_patch(circle)
    axs.axis("off")

    return sb_position, sb_size


def reconstruct_manual_general_withref_fromGUI(im, imref, sb_position, sb_size, imagesize):
    """
    

    Parameters
    ----------
    im : Object hologram
    imref : Reference hologram
    sb_position : Position of the sideband, which is found from the aforementioned function
    sb_size : Size of the applied mask, again WPO is used and 1/3 the distance from sb to cb is used
    imagesize : The size of the output reconstructed hologram. In our case due to stitching, each hologram is 
    reduced in shape by 1/16, i.e. a 4k image becomes a 256px image. 

    Returns
    -------
    im_reconstruct_P : The reconstructed phase image, which is later unwrapped.

    """
    im_reconstruct_P = im.reconstruct_phase(imref,
                                            sb_position=sb_position, sb_size=sb_size, output_shape=(imagesize[0]*1/16, imagesize[1]*1/16))
    return im_reconstruct_P


# %%
# load files
path = r'C:\Users\User\HoloData'
    
files = glob.glob(path+"\\**.dm4")
#Sort files in terms of the string, e.g. natsort sorts according to a common trend in the file string
#e.g. file_1, file_3,file_4, file_2 --> file_1,file_2,file_3,file_4,
files_sorted = natsort.natsorted(files)

#Define the group size of the acquired data
N_Size = 3
# Split File into equal size chunks depending on the acquired data. In our case, N=3 (+2V, 0V, -2V)
files_split = list(zip(*[iter(files_sorted)] * N_Size))

# save_path_unwr = r'\\ait-pdfs.win.dtu.dk\Services\NLAB\cen-archive\P87682-MolecularWindows\Mads_S_Larsen\ETEM\2023\03112023_MathiasAda_ESPPs\Holos\SavePath_1Quarter\Unwrapped'
plt.close("all")

#define empty vectors to fill
AllImgs = []
col = []
col0 = []
col_2 = []
col_3 = []
meancol_vec = []
meanrow_vec = []
#The stitching parameters with offsets. The size is predefined by knowing the final shape of the stitched image
meancol_mat = np.zeros((8,9))


count3 = 0
#Overview of counters:
    # count1 : fills up in the row, is reset to 0 whenever a new column is applied
    # count2 : fills up for the columns, is set to j-1 for the j'th column
    # count3 : is used to save count1 when the column is completed - 
    # is checked against to ensure the j-1 col is completed and then saved in array
for i, lst in enumerate(files_split): #loop through all lists defined in line 133
    if i >= 0:
        for k, item in enumerate(lst):
            if k == 0:
                data_obj = hs.load(item, signal_type="hologram")
                it0 = item.split("\\")[-1].split(".")[0].split("68kX")[-1]
            if k == 1:
                data_ref = hs.load(item, signal_type="hologram")
                itref = item.split("\\")[-1].split(".")[0].split("68kX")[-1]
            # if k==2:
            #     data_neg2VBias = hs.load(item,signal_type="hologram")
            #     it2 = item.split("\\")[-1].split(".")[0].split("68kX")[-1]
    image_size = np.shape(data_obj)
    # nm #pos y
    Stage_Y = abs(data_obj.metadata.Acquisition_instrument.TEM.Stage.y)
    scale = data_obj.axes_manager["x"].scale
    new_scale = (1/16)**-1*scale

    # nm #pos x
    Stage_X = abs(data_obj.metadata.Acquisition_instrument.TEM.Stage.x)
    if i == 0:
        Stage_Y_0 = Stage_Y  # initial y, row
        Stage_X_0 = Stage_X
        sb_poz, sbsiz = Get_sb_and_Size(data_obj)

    diff_y = Stage_Y - Stage_Y_0

    diff_x = abs((Stage_X - Stage_X_0)*1E3)

    if diff_y <= 1E-4:

        Unwrapped_Phase_vec = []  # reset for every new column
        Reconstructed_Hologram_vec = []
    Reconstructed_Hologram = reconstruct_manual_general_withref_fromGUI(
        data_obj, data_ref, sb_poz, sbsiz, image_size)

    Reconstructed_Hologram_phase = Reconstructed_Hologram.phase()
    Reconstructed_Hologram_phase = Reconstructed_Hologram_phase[10:-10, :]
    
    Unwrapped_Phase = unwrap_phase(Reconstructed_Hologram_phase)
    if i == 0:
        ImgSize = Unwrapped_Phase.shape # Size of the hologram
        black_img = np.zeros((ImgSize[0], ImgSize[1]))  # create black canvas
        
    if np.mean(Unwrapped_Phase) < 0:
        Unwrapped_Phase = Unwrapped_Phase*-1
        
    if diff_x <= 0.65:
        count2 = 0
        if diff_y <= 1E-4:
            meanrow = abs(np.mean(Unwrapped_Phase[-30:, :]))
            Unwrapped_Phase_vec += [Unwrapped_Phase]
            count1 = 0

        if diff_y >= 1E-4:
            Unwrapped_Phase = Unwrapped_Phase + meanrow
            meanrow = abs(np.mean(Unwrapped_Phase[-30:, :]))
            Unwrapped_Phase_vec += [Unwrapped_Phase]
            col0 = np.concatenate([img for img in Unwrapped_Phase_vec])
            count1 += 1
        if float(count1)==meancol_mat.shape[0]: #skip last image 
            continue
        meancol_i = abs(np.mean(Unwrapped_Phase[:, -30:]))
        meancol_mat[count1,count2] = meancol_i #.append(meancol_i)
    
    if diff_x >= 0.65:  # movement in col
        
        
        NumOfImgs = len(col0)/Unwrapped_Phase.shape[0]
        
        if diff_y < 1E-4:
            count1 = 0
            
            meancol = meancol_mat[count1,count2]
            
            count2 += 1 #add, i.e. new column added
            col = []  # reset col
            Unwrapped_Phase = Unwrapped_Phase + meancol
            Unwrapped_Phase_vec += [Unwrapped_Phase]

            meanrow = abs(np.mean(Unwrapped_Phase[-30:, :]))
            

        if diff_y >= 1E-4:  # moving downwards
            meancol = meancol_mat[count1,count2] #row,column
            
            
            Unwrapped_Phase = Unwrapped_Phase + meanrow
            
            Unwrapped_Phase = Unwrapped_Phase + meancol
            
            meanrow = abs(np.mean(Unwrapped_Phase[-30:, :]))
            
            meancol_i = abs(np.mean(Unwrapped_Phase[:, -30:]))
            
            meancol_mat[count1,count2] = meancol_i
            
            count1 += 1
            
            Unwrapped_Phase_vec += [Unwrapped_Phase]
            col = np.concatenate([img for img in Unwrapped_Phase_vec])
            print("counter = " + str(count1))
            if count1 > 0:
                diffImgs = int(NumOfImgs - count1 - 1)
                repblack = np.tile(black_img.T, diffImgs).T
                col_3 = np.concatenate([col, repblack], axis=0)

        if count1 > 0:
            count3 = count1  # set counter to a new saved counter
        if diff_y <= 1E-4:  # check for a new column
            # check and verify counter 3 is not zero (i.e., the column is done)
            if count3 > 0:
                AllImgs.append(col_3)
# %%
#Construct and hsow image
plt.close("all")
# AllImgs = np.asarray(AllImgs, dtype=float)
# AllImgs_2 = np.concatenate(AllImgs, axis=0)
AllImgs_Newer = np.zeros((8, 1888, 256))
AllImgs_Newer[0, :, :] = col0
AllImgs_Newer[1:, :, :] = AllImgs

AllImgsNewer_Offset = []
AllImgsNewer_Offset.append(col0)
for i in range(AllImgs_Newer.shape[0]-1):
    col_i = AllImgs_Newer[i]
    col_i1 = AllImgs_Newer[i+1]

    meancol = np.mean(col_i[:, -30:])
    col_i1 = col_i1 + meancol
    # print(meancol)
    AllImgsNewer_Offset.append(col_i1)

for col in AllImgsNewer_Offset:
    meancol = np.mean(col[:, -30:])
    print(meancol)
    
# plt.close("all")
Image_Final_radial = np.concatenate(AllImgs_Newer, axis=1)
# Image_Final_rec = np.concatenate([col0_rec,Img_Cols_rec],axis=1)
# Add some blur
# Image_Final = ndimage.median_filter(Image_Final,size=5)
# Image_Final = ndimage.gaussian_filter(Image_Final,sigma=3)
# Show image
Image_Final_radial = np.rot90(Image_Final_radial)
new_scale2 = new_scale  # *NumOfImgs
figunwr, axsunwr = plt.subplots(tight_layout=True)
imzhow = axsunwr.imshow(Image_Final_radial, cmap="inferno",vmin=0,vmax=30)
axsunwr.axis("off")

cbar = figunwr.colorbar(imzhow, orientation='vertical')
cbar.ax.tick_params(labelsize=14)
cbar.ax.set_ylabel('$\Delta \phi$ [rad]', fontsize=25)
font_scalebar = {'family': "Arial",
                 'size': 30,
                 'weight': 'bold'}
scalebar = ScaleBar(new_scale2,
                    "nm", length_fraction=0.5,
                    location='lower right', box_alpha=0, color='w',
                    font_properties=font_scalebar
                    )
axsunwr.add_artist(scalebar)
#%%
#Save image and arrays

save_path = path = r'C:\Users\User\HoloData\SavePath'
figunwr.savefig(
    save_path+"\\"+"FirstQuadrant_UnwrappedPhase_2V_WithOffSet_3.png")

np.save(save_path+"\\"+"FirstQuadrant_UnwrappedPhase_2V_WithOffSet_NumpyArray.npy",Image_Final_radial)
np.savetxt(save_path+"\\"+"FirstQuadrant_UnwrappedPhase_2V_WithOffSet_TXT.txt",Image_Final_radial)
