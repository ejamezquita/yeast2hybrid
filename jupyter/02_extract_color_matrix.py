import matplotlib.pyplot as plt
import numpy as np
from glob import glob
import os
from scipy import ndimage
import utils
import argparse


def main():
    parser = argparse.ArgumentParser(description='Extract a color matrix from a plate')
    parser.add_argument('src', metavar='raw_dir', type=str, help='directory where raw images are located')
    parser.add_argument('name', metavar='name_dir', type=str, help='name of plate origin')
    parser.add_argument('dst', metavar='dst_dir', type=str, help='directory where processed data are located')
    parser.add_argument('diagdst', metavar='diagdst_dir', type=str, help='directory where diagnostic images are stored')
    parser.add_argument('plate_num', metavar='plate_num', type=int, help='Number of plates to process')
    parser.add_argument('color_check', metavar='color_check', type=int, help='Color chanel to check to rotate')
    parser.add_argument('check_rotation', metavar='check_rotation', type=int, help='Check if image is rotated before loading')
    parser.add_argument('verbose', metavar='generate_diagnostic', type=int, help='Whether diagnostic images are produced')
    args = parser.parse_args()


    src = '..' + os.sep + args.src + os.sep
    lsrc = src + args.name + os.sep
    dst = '..' + os.sep + args.dst + os.sep

    genes = os.listdir(lsrc)
    size = 725
    R = 25

    for gidx in range(len(genes)):
        print('Working with gene', genes[gidx])
        gdst = dst + lsrc.split(os.sep)[-2] + os.sep + genes[gidx] + os.sep
        
        for platenum in range(len(args.plate_num)):

            platefile = glob(lsrc + genes[gidx] + os.sep + '*_{:02d}*'.format(platenum) )[0]
            bname = os.path.split(os.path.splitext(platefile)[0])[1]
            rgb = utils.load_image(platefile, color_check=1, check_rotation=False)

            filename = gdst + bname + '_plateslice.csv'
            meta = np.loadtxt(filename, delimiter=',', dtype=int)
            filename = gdst + bname + '_centers.npy'

            plateslice = np.s_[ meta[0]:meta[1], meta[2]:meta[3] ]
            diagnostic = rgb[plateslice]

            coords = np.load(filename, allow_pickle=True)
            nrows = coords.shape[0]
            ncols = coords.shape[1]

            maximg = np.max(diagnostic, axis=2)
            stdimg = np.std(diagnostic, axis=2)
            background = (maximg < meta[4]) & (stdimg < meta[5])

            img = np.where(diagnostic[:,:,1] > diagnostic[:,:,2], diagnostic[:,:,1], diagnostic[:,:,2])

            eimg = ndimage.grey_erosion(diagnostic[:,:,1], size=size, mode='constant', cval=255)
            cimg = ndimage.grey_dilation(eimg, size=size, mode='constant', cval=0)
            unif = ndimage.uniform_filter(cimg, size=size//3, mode='reflect')
            unif[ unif > img ] = img[unif > img]
            corr2 = img - unif
            corr2[background] = 0

            q = np.zeros((nrows, ncols))

            for row in range(nrows):
                for col in range(ncols):
                    pos = coords[row,col,:].astype(int)
                    if np.sum(pos) > 10:
                        slic = np.s_[pos[1] - R : pos[1] + R + 1, pos[0] - R : pos[0] + R + 1]
                        foo = corr2[slic]
                        q[row,col] = np.quantile(foo[foo > 0], 0.5)
                    else:
                        q[row,col] = 0

            filename = gdst + bname + '_colormatrix.csv'
            print(filename)
            np.savetxt(filename, q, delimiter=',', fmt='%.1f')

if args.verbose:

fig, ax = plt.subplots(1,1, figsize=(13,9), sharex=True, sharey=True)
ax = np.atleast_1d(ax).ravel()
ax[0].imshow(q, cmap='Reds_r', vmin=0)
ax[0].set_xticks(range(ncols), range(ncols))
ax[0].set_yticks(range(nrows), range(nrows));


# In[27]:


fig, ax = plt.subplots(2,2, figsize=(14,10), sharex=True, sharey=True)
ax = np.atleast_1d(ax).ravel(); i = 0

vmax = np.max(corr2)

ax[i].imshow(img, vmin=0, origin='upper', cmap='Reds_r'); i+=1
ax[i].imshow(cimg, vmin=0, vmax=255, origin='upper', cmap='inferno'); i+=1
ax[i].imshow(unif, vmin=0, vmax=255, origin='upper', cmap='inferno'); i+=1
ax[i].imshow(corr2, vmin=0, vmax=vmax, origin='upper', cmap='Reds_r'); i+=1

for i in range(len(ax)):
    ax[i].axis('off')
    
fig.tight_layout();


# In[109]:


row = 11
fig, ax = plt.subplots(6,8, figsize=(20,15), sharex=True, sharey=True)
ax = np.atleast_1d(ax).ravel()

for i in range(len(ax)):
    pos = coords[row,i]
    if np.sum(pos) > 10:
        slic = np.s_[pos[1] - R : pos[1] + R + 1, pos[0] - R : pos[0] + R + 1]
        foo = corr2[slic]
        med = np.median(foo[foo > 0])
        ax[i].imshow(foo, origin='upper', cmap='Reds_r', vmax=vmax, vmin=0);
        ax[i].set_title(int(med))
    ax[i].axis('off')

fig.tight_layout();

if __name__ == '__main__':
    main()
