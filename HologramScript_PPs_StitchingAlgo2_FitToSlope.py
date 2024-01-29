# -*- coding: utf-8 -*-
"""
Created on Fri Jan 26 13:54:02 2024

@author: maslar
"""

import matplotlib.pyplot as plt
import numpy as np
import glob
from mpl_toolkits.axes_grid1 import make_axes_locatable

scale = 2.1154134273529053

#%%
path = r'C:\Users\User\HoloFiles\Temporary_SavePath'
file = glob.glob(path+"\\**.npy")[0]
plt.close("all")
# AllImgs = np.asarray(AllImgs, dtype=float)
# AllImgs_2 = np.concatenate(AllImgs, axis=0)


plt.close("all")    
# plt.close("all")
Image_Final_radial = np.load(file)
fig,axs = plt.subplots()
axs.imshow(Image_Final_radial)
ImgSize = (236,256)
slopes_x=[]
slopes_y = []
for i in range(int(Image_Final_radial.shape[0]/ImgSize[0])):
    for j in range(int(Image_Final_radial.shape[1]/ImgSize[1])):
        Single_hol = Image_Final_radial[i*ImgSize[0]:(i+1)*ImgSize[0],j*ImgSize[1]:(j+1)*ImgSize[1]]
        mean_Single_hol = np.mean(Single_hol)
        
        Single_hol = Single_hol - mean_Single_hol
        slope_x = np.polyfit(np.linspace(0,Single_hol.shape[1],Single_hol.shape[1]),np.mean(Single_hol,axis=0),1)
        slope_y = np.polyfit(np.linspace(0,Single_hol.shape[0],Single_hol.shape[0]),np.mean(Single_hol,axis=1),1)
        
        if slope_x[0] == 0:
            continue
        slopes_x.append(slope_x)
        slopes_y.append(slope_y)
slopes_x = np.asarray(slopes_x)
slopes_y = np.asarray(slopes_y)
slopes_x_m = np.mean(slopes_x[:,0])
slopes_y_m = np.mean(slopes_y[:,0])

Img_Slope_Fixed = np.zeros((Image_Final_radial.shape[0],Image_Final_radial.shape[1]))
Img_Zero = np.zeros((Image_Final_radial.shape[0],Image_Final_radial.shape[1]))
Slope_Img_All = np.zeros((Image_Final_radial.shape[0],Image_Final_radial.shape[1]))
Slope_img_y_v = []
for j in range(int(Image_Final_radial.shape[1]/ImgSize[1])):
    x = np.linspace(j*ImgSize[1], (j+1)*ImgSize[1]-1, ImgSize[1])*slopes_x_m
    Slope_img_x = np.tile(x, (ImgSize[0], 1))

    
    for i in range(int(Image_Final_radial.shape[0]/ImgSize[0])):
        
        
        y = np.linspace(i*ImgSize[0],(i+1)*ImgSize[0]-1,ImgSize[0])*slopes_y_m
        Slope_img_y = np.tile(y,(ImgSize[1],1)).T
        
        Single_hol = Image_Final_radial[i*ImgSize[0]:(i+1)*ImgSize[0],j*ImgSize[1]:(j+1)*ImgSize[1]]
        
        mean_Single_hol = np.mean(Single_hol)
        
        Single_hol_sub = Single_hol - mean_Single_hol
        
        Img_Zero[i*ImgSize[0]:(i+1)*ImgSize[0],j*ImgSize[1]:(j+1)*ImgSize[1]] = Single_hol_sub
        
        Slope_Img_fin = Slope_img_x + Slope_img_y
        
        
        if mean_Single_hol != 0:
            Slope_Img_All[i*ImgSize[0]:(i+1)*ImgSize[0],j*ImgSize[1]:(j+1)*ImgSize[1]] = Slope_Img_fin
            Single_hol_slopeadd = Single_hol_sub + Slope_Img_fin#Slope_img + Slope_img2 
        elif mean_Single_hol == 0:
            Slope_Img_All[i*ImgSize[0]:(i+1)*ImgSize[0],j*ImgSize[1]:(j+1)*ImgSize[1]] = Single_hol
            Single_hol_slopeadd = Single_hol
            
        Img_Slope_Fixed[i*ImgSize[0]:(i+1)*ImgSize[0],j*ImgSize[1]:(j+1)*ImgSize[1]] = Single_hol_slopeadd
        
        offset_slope_y = np.mean(Slope_img_y[:,:],axis=1)[-1]
        Slope_img_y = Slope_img_y - np.mean(Slope_img_y,axis=1)[0] + offset_slope_y
    offset_slope_x = np.mean(Slope_img_x,axis=0)[-1]
    Slope_img_x = Slope_img_x - np.mean(Slope_img_x,axis=0)[0] + offset_slope_x

fig, axs = plt.subplots(tight_layout=True,figsize=(6*1.33*2,6),ncols=2)
imzhowz = []
imzhowz.append(axs[0].imshow(Img_Slope_Fixed, cmap="inferno"))#,vmin=-5,vmax=np.max(Img_Slope_Fixed))
axs[0].set_title("2D Phase Image",fontsize=20)


imzhowz.append(axs[1].imshow(Slope_Img_All, cmap="inferno"))#,vmin=-5,vmax=np.max(Img_Slope_Fixed))
axs[1].set_title("2D Slope",fontsize=20)

for n in range(2):
    axs[n].axis("off")
    divider = make_axes_locatable(axs[n])
    cax = divider.append_axes("bottom", size="5%", pad=0.25)
     
    # Plot horizontal colorbar
    cbar = plt.colorbar(imzhowz[n], orientation="horizontal", cax=cax)
    cbar.ax.tick_params(labelsize=18)
    cbar.ax.set_xlabel('$\Delta \phi$ [rad]', fontsize=25)

    axs[n].tick_params(axis="both",labelsize=16)
    
save_path = r'C:\Users\User\HoloFiles\SavePath'
fig.savefig(
    save_path+"\\"+"FirstQuadrant_UnwrappedPhase.png",bbox_inches="tight")
np.save(save_path+"\\"+"FirstQuadrant_UnwrappedPhas_Numpy.npy",Img_Slope_Fixed)
np.savetxt(save_path+"\\"+"FirstQuadrant_UnwrappedPhase_TXT.txt",Img_Slope_Fixed)

#%%
plt.close("all")
path = r'C:\Users\maslar\OneDrive - Danmarks Tekniske Universitet\MowinShared\Project PhD Mads\Mathias_Ada_PhasePlates\Temporary_Save\Linescans'
save_path = r'C:\Users\maslar\OneDrive - Danmarks Tekniske Universitet\MowinShared\Project PhD Mads\Mathias_Ada_PhasePlates\New_Approach_Slopes'
files = glob.glob(path+"\\**.npy")
offsets_elec = [0]
mean_vec_all = []
x_vec_all = []
for i,file in enumerate(files):
    Image = np.load(file)
    slopes_y = []
    for j in range(int(Image.shape[0]/ImgSize[0])):
        
        Single_hol = Image[j*ImgSize[0]:(j+1)*ImgSize[0],:]
        
        mean_Single_hol = np.mean(Single_hol)
        
        Single_hol = Single_hol - mean_Single_hol
        
        
        slope_y = np.polyfit(np.linspace(0,Single_hol.shape[0],Single_hol.shape[0]),np.mean(Single_hol,axis=1),1)
        if slope_y[0] < 0:
            Single_hol = Single_hol*-1
        slope_y = np.polyfit(np.linspace(0,Single_hol.shape[0],Single_hol.shape[0]),np.mean(Single_hol,axis=1),1)
        slopes_y.append(slope_y)
    slopes_y = np.asarray(slopes_y)
    mean_slope_y = np.mean(slopes_y[:,0])
    mean_slope_y = mean_slope_y 
    fig,axs = plt.subplots(tight_layout=True)
    Slope_Imgs = []
    mean_vec = []
    x_vec = []
    for j in range(int(Image.shape[0]/ImgSize[0])):
        y = np.linspace(j*ImgSize[0],(j+1)*ImgSize[0]-1,ImgSize[0])*mean_slope_y + offsets_elec[i]
        Slope_img_y = np.tile(y,(ImgSize[1],1)).T 
        
        Slope_Imgs.append(Slope_img_y)
        Single_hol = Image[j*ImgSize[0]:(j+1)*ImgSize[0],:] 
         
        mean_Single_hol = np.mean(Single_hol)
        Single_hol = Single_hol - mean_Single_hol
        
        slope_y = np.polyfit(np.linspace(0,Single_hol.shape[0],Single_hol.shape[0]),np.mean(Single_hol,axis=1),1)
        if slope_y[0] < 0:
            Single_hol = Single_hol*-1

        hologram_w_slope = Single_hol  + Slope_img_y 
        
        mean = np.mean(hologram_w_slope, axis=1)   
        x = np.linspace(0,len(mean),len(mean))
        
        
        if j==int(Image.shape[0]/ImgSize[0]-1):
            offsets_elec.append(mean[-1])
        mean_vec.append(mean)
        x_vec.append(x)
        axs.plot(x+j*x[-1],mean)
        axs.plot(x+j*x[-1],np.polyval(np.polyfit(x,mean,1),x),'k--')
    Slope_Imgs_2 = np.concatenate([img for img in Slope_Imgs], axis=0)
    
    fig,axs = plt.subplots()
    axs.imshow(Slope_Imgs_2)
    mean_vec_2 = np.concatenate([meang for meang in mean_vec])
    x_vec_2 = np.concatenate([xi for xi in x_vec])
    mean_vec_all.append(mean_vec_2)
    x_vec_all.append(x_vec_2)
    file_name = file.split("\\")[-1].split("NotFixed")[0]
    print(file_name)
    np.save(save_path+"\\"+file_name+"_NumpyArray.npy",mean_vec_2)
    np.savetxt(save_path+"\\"+file_name+"_TXT.txt",mean_vec_2)
    
    