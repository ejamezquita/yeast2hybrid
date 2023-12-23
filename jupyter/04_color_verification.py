import matplotlib.pyplot as plt

import numpy as np
from glob import glob
import os

import pandas as pd
import utils

import argparse

"""
How to run and expected output

$ python3 04_color_verification.py raw proc leyre mather p53 diagnostic 6
Reference:	p53
Colonies to plot:	['LBD37']
Comparing gene LBD37 with reference p53
Will save all results in  ../diagnostic/leyre/LBD37/
Generate image ../diagnostic/leyre/LBD37/p53_vs_20231130_plate_01_sigma_significant.jpg
Generate image ../diagnostic/leyre/LBD37/p53_vs_20231130_plate_01_laplace_significant.jpg
Generate image ../diagnostic/leyre/LBD37/p53_vs_20231130_plate_02_sigma_significant.jpg
Generate image ../diagnostic/leyre/LBD37/p53_vs_20231130_plate_02_laplace_significant.jpg
Generate image ../diagnostic/leyre/LBD37/p53_vs_20231130_plate_03_sigma_significant.jpg
Generate image ../diagnostic/leyre/LBD37/p53_vs_20231130_plate_03_laplace_significant.jpg
Generate image ../diagnostic/leyre/LBD37/p53_vs_20231130_plate_04_sigma_significant.jpg
Generate image ../diagnostic/leyre/LBD37/p53_vs_20231130_plate_04_laplace_significant.jpg
Generate image ../diagnostic/leyre/LBD37/p53_vs_20231130_plate_05_sigma_significant.jpg
Generate image ../diagnostic/leyre/LBD37/p53_vs_20231130_plate_05_laplace_significant.jpg
Generate image ../diagnostic/leyre/LBD37/p53_vs_20231130_plate_06_sigma_significant.jpg
Generate image ../diagnostic/leyre/LBD37/p53_vs_20231130_plate_06_laplace_significant.jpg
"""

def main():
    parser = argparse.ArgumentParser(description='Extract a color matrix from a plate')
    parser.add_argument('rsrc', metavar='raw_dir', type=str, help='directory where raw images are located')
    parser.add_argument('psrc', metavar='raw_dir', type=str, help='directory where processed images are located')
    parser.add_argument('name', metavar='name_dir', type=str, help='name of plate origin')
    parser.add_argument('refdir', metavar='refname_dir', type=str, help='name of reference subfolder')
    parser.add_argument('refgene', metavar='refname_dir', type=str, help='name of reference gene')
    parser.add_argument('diagdst', metavar='diagdst_dir', type=str, help='directory where diagnostic images are stored')
    parser.add_argument('plate_num', metavar='plate_num', type=int, help='Number of plates to process')
    args = parser.parse_args()

    reference_foldername = args.refdir
    reference_gene = args.refgene

    fs = 15
    R = 80
    figR = 2
    axcolmax = 10

    src = '..' + os.sep + args.psrc + os.sep
    lsrc = src + args.name + os.sep
    msrc = src + reference_foldername + os.sep + reference_gene + os.sep

    rsrc = '..' + os.sep + args.rsrc + os.sep 

    dst = '..' + os.sep + args.diagdst + os.sep

    genes = os.listdir(lsrc)
    print('Reference:\t', reference_gene, '\nColonies to plot:\t', genes, sep='')


    for gidx in range(len(genes)):
        gsrc = lsrc + genes[gidx] + os.sep
        ddst = dst + lsrc.split(os.sep)[-2] + os.sep + genes[gidx] + os.sep

        print('Comparing gene', genes[gidx], 'with reference', reference_gene )
        print('Will save all results in ',ddst)

        for platenum in range(1, args.plate_num +1):
            sigma_file = glob(gsrc + reference_gene + '*_{:02d}_sigma_significant_differences.csv'.format(platenum) )[0]
            sbname = '_'.join((os.path.split(os.path.splitext(sigma_file)[0])[1]).split('_')[:-1])

            laplace_file = glob(gsrc + reference_gene + '*_{:02d}_laplace_significant_differences.csv'.format(platenum) )[0]
            lbname = '_'.join((os.path.split(os.path.splitext(laplace_file)[0])[1]).split('_')[:-1])

            ssignif = pd.read_csv(sigma_file, index_col=0)
            lsignif = pd.read_csv(laplace_file, index_col=0)

            rawimg_file = rsrc + reference_foldername + os.sep + reference_gene + os.sep
            rawimg_file = glob(rawimg_file + '*_{:02d}*'.format(platenum))[0]
            rname = os.path.split(os.path.splitext(rawimg_file)[0])[1]

            pltimg_file = gsrc.replace(src,rsrc)
            pltimg_file = glob(pltimg_file + '*_{:02d}*'.format(platenum))[0]
            bname = os.path.split(os.path.splitext(pltimg_file)[0])[1]

            filename = msrc + rname + '_plateslice.csv'
            meta = np.loadtxt(filename, delimiter=',', dtype=int)
            plateslice = np.s_[ meta[0]:meta[1], meta[2]:meta[3] ]
            rawimg = utils.load_image(rawimg_file, check_rotation=meta[14], color_check=meta[15])
            rawimg = rawimg[plateslice]
            filename = msrc + rname + '_centers.npy'
            rawcoords = np.load(filename, allow_pickle=True)

            filename = gsrc + bname + '_plateslice.csv'
            meta = np.loadtxt(filename, delimiter=',', dtype=int)
            plateslice = np.s_[ meta[0]:meta[1], meta[2]:meta[3] ]
            pltimg = utils.load_image(pltimg_file, check_rotation=meta[14], color_check=meta[15])
            pltimg = pltimg[plateslice]
            filename = gsrc + bname + '_centers.npy'
            pltcoords = np.load(filename, allow_pickle=True)

            filename = ddst + sbname + '.jpg'
            if not os.path.isfile(filename):

                axrows = min([4, len(ssignif)//axcolmax + 1])
                axcols = min([axcolmax, len(ssignif)])
                plotnums = min([axrows*axcolmax, len(ssignif)])

                fig, ax = plt.subplots(2*axrows, axcols, figsize=(figR*axcols-1, 2*figR*axrows), sharex=True, sharey=True)
                ax = ax.reshape(2*axrows, axcols)

                for i in range(plotnums):
                    ix, jx = i//axcols, i%axcols
                    ax[2*ix, jx].set_title(ssignif.iloc[i,0], fontsize=fs)
                    ax[2*ix+1, jx].set_title('$\Delta I =${:.2f}$\\sigma$'.format(ssignif.iloc[i,7]), fontsize=fs)

                    row, col = [int(foo[1:])-1 for foo in ssignif.iloc[i,1].split('-')[1:] ] 

                    j = 2*ix
                    for (cc,img) in zip([rawcoords, pltcoords],[rawimg,pltimg]):
                        dots = cc[2*row : 2*row + 2, 2*col : 2*col + 2]
                        center = np.mean(np.mean(dots, axis = 0), axis = 0).astype(int)
                        rss = np.s_[center[1] - R : center[1] + R, center[0] - R : center[0] + R]
                        ax[j,jx].imshow(img[rss], vmin=0, origin='upper'); j+=1

                for i in range(1, 1 + axrows*axcols - plotnums):
                    fig.delaxes(ax[-1, -i])
                    fig.delaxes(ax[-2, -i])

                for i in range(axrows):
                    ax[2*i,0].set_ylabel(reference_gene, fontsize=fs)
                    ax[2*i + 1,0].set_ylabel(genes[gidx], fontsize=fs)

                for a in ax.ravel():
                    a.tick_params(bottom=False, left=False, labelbottom=False, labelleft=False)

                fig.suptitle(bname + ' : Significances according to St.Dev criteria', fontsize=fs+5)
                fig.tight_layout()

                print('Generate image', filename)
                plt.savefig(filename, dpi=72, format='jpg', bbox_inches='tight', pil_kwargs={'optimize':True})
                plt.close()

            filename = ddst + lbname + '.jpg'
            if not os.path.isfile(filename):

                axrows = min([4, len(lsignif)//axcolmax + 1])
                axcols = min([axcolmax, len(lsignif)])
                plotnums = min([axrows*axcolmax, len(lsignif)])

                fig, ax = plt.subplots(2*axrows, axcols, figsize=(figR*axcols-1, 2*figR*axrows), sharex=True, sharey=True)
                ax = ax.reshape(2*axrows, axcols)

                for i in range(plotnums):
                    ix, jx = i//axcols, i%axcols
                    ax[2*ix, jx].set_title(lsignif.iloc[i,0], fontsize=fs)
                    ax[2*ix+1, jx].set_title('$\\alpha(\\Delta I) =${:.3f}'.format(lsignif.iloc[i,8]), fontsize=fs)

                    row, col = [int(foo[1:])-1 for foo in lsignif.iloc[i,1].split('-')[1:] ] 

                    j = 2*ix
                    for (cc,img) in zip([rawcoords, pltcoords],[rawimg,pltimg]):
                        dots = cc[2*row : 2*row + 2, 2*col : 2*col + 2]
                        center = np.mean(np.mean(dots, axis = 0), axis = 0).astype(int)
                        rss = np.s_[center[1] - R : center[1] + R, center[0] - R : center[0] + R]
                        ax[j,jx].imshow(img[rss], vmin=0, origin='upper'); j+=1

                for i in range(1, 1 + axrows*axcols - plotnums):
                    fig.delaxes(ax[-1, -i])
                    fig.delaxes(ax[-2, -i])

                for i in range(axrows):
                    ax[2*i,0].set_ylabel(reference_gene, fontsize=fs)
                    ax[2*i + 1,0].set_ylabel(genes[gidx], fontsize=fs)

                for a in ax.ravel():
                    a.tick_params(bottom=False, left=False, labelbottom=False, labelleft=False)

                fig.suptitle(bname + ' : Significances according to Laplace assumption', fontsize=fs+5)
                fig.tight_layout()
                print('Generate image', filename)
                plt.savefig(filename, dpi=72, format='jpg', bbox_inches='tight', pil_kwargs={'optimize':True})
                plt.close()
    
if __name__ == '__main__':
    main()


