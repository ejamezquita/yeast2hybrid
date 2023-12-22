import matplotlib.pyplot as plt
import numpy as np
from glob import glob
import os
from scipy import ndimage
import utils
import argparse

"""
$ python3 02_extract_color_matrix.py raw mather proc diagnostic 6 1

Working with gene p53
Ready to process ../raw/mather/p53/p53_plate_01.jpg
../proc/mather/p53/p53_plate_01_colormatrix.csv
Generated image ../diagnostic/mather/p53/p53_plate_01_light_correction.jpg
Ready to process ../raw/mather/p53/p53_plate_02.jpg
../proc/mather/p53/p53_plate_02_colormatrix.csv
Generated image ../diagnostic/mather/p53/p53_plate_02_light_correction.jpg
Ready to process ../raw/mather/p53/p53_plate_03.jpg
../proc/mather/p53/p53_plate_03_colormatrix.csv
Generated image ../diagnostic/mather/p53/p53_plate_03_light_correction.jpg
Ready to process ../raw/mather/p53/p53_plate_04.jpg
../proc/mather/p53/p53_plate_04_colormatrix.csv
Generated image ../diagnostic/mather/p53/p53_plate_04_light_correction.jpg
Ready to process ../raw/mather/p53/p53_plate_05.jpg
../proc/mather/p53/p53_plate_05_colormatrix.csv
Generated image ../diagnostic/mather/p53/p53_plate_05_light_correction.jpg
Ready to process ../raw/mather/p53/p53_plate_06.jpg
../proc/mather/p53/p53_plate_06_colormatrix.csv
Generated image ../diagnostic/mather/p53/p53_plate_06_light_correction.jpg
"""

def main():
    parser = argparse.ArgumentParser(description='Extract a color matrix from a plate')
    parser.add_argument('src', metavar='raw_dir', type=str, help='directory where raw images are located')
    parser.add_argument('name', metavar='name_dir', type=str, help='name of plate origin')
    parser.add_argument('dst', metavar='dst_dir', type=str, help='directory where processed data are located')
    parser.add_argument('diagdst', metavar='diagdst_dir', type=str, help='directory where diagnostic images are stored')
    parser.add_argument('plate_num', metavar='plate_num', type=int, help='Number of plates to process')
    parser.add_argument('verbose', metavar='generate_diagnostic', type=int, help='Whether diagnostic images are produced')
    args = parser.parse_args()


    src = '..' + os.sep + args.src + os.sep
    lsrc = src + args.name + os.sep
    dst = '..' + os.sep + args.dst + os.sep
    diagdst = '..' + os.sep + args.diagdst + os.sep

    genes = os.listdir(lsrc)
    size = 725
    R = 25
    fs = 12

    for gidx in range(len(genes)):
        print('Working with gene', genes[gidx])
        gdst = dst + lsrc.split(os.sep)[-2] + os.sep + genes[gidx] + os.sep
        
        ddst = diagdst + lsrc.split(os.sep)[-2] + os.sep
        if not os.path.isdir(ddst):
            os.mkdir(ddst)
        ddst += genes[gidx] + os.sep
        if not os.path.isdir(ddst):
            os.mkdir(ddst)
        
        for platenum in range(1,args.plate_num+1):
                
            platefile = glob(lsrc + genes[gidx] + os.sep + '*_{:02d}*'.format(platenum) )[0]
            print('Ready to process', platefile)
            bname = os.path.split(os.path.splitext(platefile)[0])[1]
            filename = gdst + bname + '_plateslice.csv'
            meta = np.loadtxt(filename, delimiter=',', dtype=int)
            
            rgb = utils.load_image(platefile, check_rotation=meta[14], color_check=meta[15])
            
            if ~os.path.isfile(gdst + bname + '_colormatrix.csv'):

                plateslice = np.s_[ meta[0]:meta[1], meta[2]:meta[3] ]
                diagnostic = rgb[plateslice]
                
                filename = gdst + bname + '_centers.npy'
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

                    fig, ax = plt.subplots(2,2, figsize=(9,7), sharex=False, sharey=False)
                    ax = np.atleast_1d(ax).ravel(); i = 0

                    vmax = np.max(corr2)

                    ax[i].imshow(img, vmin=0, origin='upper', cmap='Reds_r'); 
                    ax[i].set_xlabel('Original Grayscale', fontsize=fs); i+=1
                    ax[i].imshow(unif, vmin=0, vmax=255, origin='upper', cmap='inferno');
                    ax[i].set_xlabel('Ilumination correction', fontsize=fs); i+=1
                    ax[i].imshow(corr2, vmin=0, vmax=vmax, origin='upper', cmap='Reds_r');
                    ax[i].set_xlabel('Adjusted Grayscale - Background', fontsize=fs); i+=1
                    ax[i].imshow(q, origin='upper', cmap='Reds_r');
                    ax[i].set_xlabel('Intensity matrix', fontsize=fs); i+=1

                    for i in range(len(ax)):
                        ax[i].tick_params(bottom=False, left=False, labelbottom=False, labelleft=False)

                    fig.suptitle(genes[gidx] + ' : ' + bname + ' : colormatrix extraction', fontsize=fs+5)
                    fig.tight_layout();
                    filename = ddst + 'light_correction_' + bname + '.jpg'
                    plt.savefig(filename, format='jpg', dpi=100, bbox_inches='tight', pil_kwargs={'optimize':True})
                    plt.close()
                    print('Generated image',filename)

if __name__ == '__main__':
    main()
