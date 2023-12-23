import matplotlib.pyplot as plt

import numpy as np
from glob import glob
import os

import pandas as pd
from importlib import reload
from skimage import exposure

import argparse

"""
How to run and expected output

$ python3 03_compare_color_matrices.py proc leyre proc mather p53 diagnostic 6 1

Reference:	p53
Genes to test:	['LBD37']
Comparing gene LBD37 with reference p53
Will save all results in  ../proc/leyre/LBD37/
../proc/leyre/LBD37/p53_vs_20231130_plate_01_laplace_significant_differences.csv
../proc/leyre/LBD37/p53_vs_20231130_plate_01_sigma_significant_differences.csv
Saved image ../diagnostic/leyre/LBD37/difference_detection_summary_20231130_plate_01.jpg
../proc/leyre/LBD37/p53_vs_20231130_plate_02_laplace_significant_differences.csv
../proc/leyre/LBD37/p53_vs_20231130_plate_02_sigma_significant_differences.csv
Saved image ../diagnostic/leyre/LBD37/difference_detection_summary_20231130_plate_02.jpg
../proc/leyre/LBD37/p53_vs_20231130_plate_04_laplace_significant_differences.csv
../proc/leyre/LBD37/p53_vs_20231130_plate_04_sigma_significant_differences.csv
Saved image ../diagnostic/leyre/LBD37/difference_detection_summary_20231130_plate_04.jpg
../proc/leyre/LBD37/p53_vs_20231130_plate_05_laplace_significant_differences.csv
../proc/leyre/LBD37/p53_vs_20231130_plate_05_sigma_significant_differences.csv
Saved image ../diagnostic/leyre/LBD37/difference_detection_summary_20231130_plate_05.jpg
../proc/leyre/LBD37/p53_vs_20231130_plate_06_laplace_significant_differences.csv
../proc/leyre/LBD37/p53_vs_20231130_plate_06_sigma_significant_differences.csv
Saved image ../diagnostic/leyre/LBD37/difference_detection_summary_20231130_plate_06.jpg
"""

def main():
    parser = argparse.ArgumentParser(description='Extract a color matrix from a plate')
    parser.add_argument('src', metavar='raw_dir', type=str, help='directory where processed images are located')
    parser.add_argument('name', metavar='name_dir', type=str, help='name of plate origin')
    parser.add_argument('dst', metavar='dst_dir', type=str, help='directory where processed data are located')
    parser.add_argument('refdir', metavar='refname_dir', type=str, help='name of reference subfolder')
    parser.add_argument('refgene', metavar='refname_dir', type=str, help='name of reference gene')
    parser.add_argument('diagdst', metavar='diagdst_dir', type=str, help='directory where diagnostic images are stored')
    parser.add_argument('plate_num', metavar='plate_num', type=int, help='Number of plates to process')
    parser.add_argument('verbose', metavar='generate_diagnostic', type=int, help='Whether diagnostic images are produced')
    args = parser.parse_args()


    fs = 15

    reference_foldername = args.refdir
    reference_gene = args.refgene
    data = pd.read_csv('..' + os.sep + 'library_eY2H.csv')

    src = '..' + os.sep + args.src + os.sep
    lsrc = src + args.name + os.sep

    msrc = src + reference_foldername + os.sep + reference_gene + os.sep
    dst = '..' + os.sep + args.dst + os.sep
    diagdst = '..' + os.sep + args.diagdst + os.sep

    genes = os.listdir(lsrc)
    print('Reference:\t', reference_gene, '\nGenes to test:\t', genes, sep='')

    for gidx in range(len(genes)):
        print('Comparing gene', genes[gidx], 'with reference', reference_gene )
        gdst = dst + lsrc.split(os.sep)[-2] + os.sep + genes[gidx] + os.sep
        print('Will save all results in ',gdst)
        ddst = diagdst + lsrc.split(os.sep)[-2] + os.sep + genes[gidx] + os.sep


    for platenum in range(1, args.plate_num+1):
        ref_file = glob(msrc + '*_{:02d}_colormatrix.csv'.format(platenum) )[0]
        plate_file = glob(lsrc + genes[gidx] + os.sep + '*_{:02d}_colormatrix.csv'.format(platenum) )[0]
        bname = '_'.join((os.path.split(os.path.splitext(plate_file)[0])[1]).split('_')[:-1])
        
        filename = gdst + reference_gene + '_vs_' + bname + '_laplace_significant_differences.csv'
        if not os.path.isfile(filename):

            reference = np.loadtxt(ref_file, delimiter=',')
            plate = np.loadtxt(plate_file, delimiter=',')
            nrows, ncols = plate.shape
            nonzeros = plate != 0

            reference[~nonzeros] = 0
            matched = exposure.match_histograms(plate, reference, channel_axis=None)

            diff = reference - matched
            diff[~nonzeros] = 0
            vlim = np.max( np.abs( [ np.min(diff), np.max(diff) ] ) )

            rhist, bins = np.histogram(reference, bins=range(1,257))
            phist, _ = np.histogram(plate, bins=bins)
            mhist, _ = np.histogram(matched, bins=bins)

            imgsize = np.sum(nonzeros)
            rcumsum = np.cumsum(rhist)/reference.size
            pcumsum = np.cumsum(phist)/imgsize
            mcumsum = np.cumsum(mhist)/imgsize

            means = np.zeros((nrows//2, ncols//2))
            for i in range(0,nrows,2):
                for j in range(0, ncols, 2):
                    means[i//2,j//2] = np.mean(diff[i:i+2, j:j+2])

            dhist, dbins = np.histogram(diff[nonzeros], bins=np.linspace(-vlim-1, vlim+1, 101), density=True)
            dcumsum = np.cumsum(dhist)
            dcdf = (dcumsum - dcumsum[0])/(dcumsum[-1] - dcumsum[0])
            xvals = dbins[:-1].copy()

            meanhist, _ = np.histogram(means, bins=dbins, density=True)
            meancumsum = np.cumsum(meanhist)
            meancdf = (meancumsum - meancumsum[0])/(meancumsum[-1] - meancumsum[0])

            mu, sigma = np.mean(diff), np.std(diff)
            normal = 1/(np.std(diff) * np.sqrt(2 * np.pi))*np.exp( - (xvals - np.mean(diff))**2 / (2 * np.std(diff)**2))
            ncumsum = np.cumsum(normal)
            ncdf = (ncumsum - ncumsum[0])/(ncumsum[-1] - ncumsum[0])

            b = np.std(diff)/np.sqrt(2)
            laplace = np.exp(-abs(xvals)/b)/(2.*b)

            lcdf = np.zeros(len(xvals))
            lcdf[xvals < 0] = 0.5*np.exp((xvals[xvals < 0])/b)
            lcdf[xvals >= 0] = 1-0.5*np.exp(-(xvals[xvals >= 0])/b)

            # FIND OUTLIERS BY SOME CRITERIA

            # 1. Assume differences follow a Laplace distrubution. Outliers will be those outside the 1-alpha interval
            alpha = 0.025
            lthreshold = np.abs(b * np.log(2*alpha))
            ltmask = np.abs(diff) > lthreshold

            # 2. Outliers are those 3 stds away from the mean
            sigmalim=2
            sthreshold = sigmalim*np.std(diff)
            stmask = np.abs(diff) > sthreshold

            def get_significant_summary(tmask):
                signifscores = np.zeros((nrows//2, ncols//2), dtype=np.uint8)

                for i in range(0,nrows,2):
                    for j in range(0, ncols, 2):
                        signifscores[i//2, j//2] = np.sum(tmask[i:i+2, j:j+2])
                
                signifcoords = np.asarray(np.nonzero(signifscores > 2))
                
                signifidx = np.zeros(len(signifcoords[0]), dtype=int) - 1
                diffvals = np.zeros(len(signifidx))
                
                for i in range(len(signifidx)):
                    coordinate = 'p{:02d}-r{:02d}-c{:02d}'.format(platenum, signifcoords[0,i]+1, signifcoords[1,i]+1)
                    ref = data[data['Coordinate'] == coordinate]
                    if len(ref) > 0:
                        signifidx[i] = ref.index[0]
                        diffvals[i] = means[signifcoords[0,i], signifcoords[1,i]]
                
                signifmask = signifidx > -1
                
                signif = data.loc[signifidx[signifmask]]
                signif['IntensityDiff'] = diffvals[signifmask]
                signif['SigmaValue'] = diffvals[signifmask]/sigma
                signif['AlphaValue'] = np.where(diffvals[signifmask] < 0, .5*np.exp(diffvals[signifmask]/b), 1-.5*np.exp(-diffvals[signifmask]/b) )
                signif['ReferenceColor'] = np.where(diffvals[signifmask] < 0, 'redish', 'whiteish')
                signif['ControlColor'] = np.where(diffvals[signifmask] < 0, 'whiteish', 'redish')
                signif.sort_values(by='IntensityDiff', axis=0, ascending=False, key=np.abs, inplace=True)

                return signif


            lsignif = get_significant_summary(ltmask)
            ssignif = get_significant_summary(stmask)

            filename = gdst + reference_gene + '_vs_' + bname + '_laplace_significant_differences.csv'
            print(filename)
            lsignif.to_csv(filename, index=True, index_label='OriginalIndex')

            filename = gdst + reference_gene + '_vs_' + bname + '_sigma_significant_differences.csv'
            print(filename)
            ssignif.to_csv(filename, index=True, index_label='OriginalIndex')
            
            if args.verbose:

                fig, ax = plt.subplots(2,3, figsize=(12,7), sharex=False, sharey=False)
                ax = np.atleast_1d(ax).ravel(); i = 0

                ax[i].plot(bins[:-1], rhist, c='blue', ds='steps', label='reference')
                ax[i].plot(bins[:-1], phist, c='red', ds='steps', label='plate')
                ax[i].plot(bins[:-1], mhist, c='gray', ds='steps', alpha=0.5, label='matched')
                ax[i].set_title('Intensity distribution', fontsize=fs)

                i+=1
                ax[i].margins(0)
                ax[i].plot(xvals, dhist, c='k', lw=2, ds='steps-post', label='differences', zorder=1)
                ax[i].plot(xvals, normal, c='red', lw=2, label='N(0,{:.1f})'.format(sigma), alpha=0.5 , zorder = 3)
                ax[i].plot(xvals, laplace, c='blue', lw=2, label='L(0,{:.1f})'.format(b), alpha=0.5, zorder = 4 )
                ax[i].set_title('Probability Distribution Function', fontsize=fs)

                i+=1
                ax[i].imshow(diff*ltmask, cmap='seismic', vmin=-vlim, vmax=vlim)
                ax[i].set_xlabel('Keeping {} groups'.format(len(lsignif)), fontsize=fs)
                ax[i].set_title('At $\\alpha=$ {}, discard $|\\Delta I| <$ {:.1f}'.format(alpha, lthreshold), fontsize=fs)

                i+=1
                ax[i].plot(bins[:-1], rcumsum, c='blue', ds='steps', label='reference')
                ax[i].plot(bins[:-1], pcumsum, c='red', ds='steps', label='plate')
                ax[i].plot(bins[:-1], mcumsum, c='k', ds='steps', label='matched')
                ax[i].set_title('Cumulative sum of intensities', fontsize=fs)

                i+=1
                ax[i].plot(xvals, dcdf, c='k', label='differences', ds='steps', lw=2, alpha=1, zorder=1)
                ax[i].plot(xvals, lcdf, c='r', label='N(0,{:.1f})'.format(sigma), lw=2, alpha=0.5, zorder=3)
                ax[i].plot(xvals, ncdf, c='b', label='L(0,{:.1f})'.format(b), lw=2, alpha=0.5, zorder=4)
                ax[i].set_title('Cumulative Distribution Function', fontsize=fs)

                i+=1
                ax[i].imshow(diff*stmask, cmap='seismic', vmin=-vlim, vmax=vlim)
                ax[i].set_xlabel('Keeping {} groups'.format(len(ssignif)), fontsize=fs)
                ax[i].set_title('At $\\Delta\\sigma=$ {}, discard $|\\Delta I| <$ {:.1f}'.format(sigmalim, sthreshold), fontsize=fs)

                for i in [0,3]:
                    ax[i].set_xlabel('$I$', fontsize=fs)
                for i in [1,4]:
                    ax[i].set_xlabel('$\\Delta I$', fontsize=fs)
                for i in [2,5]:
                    ax[i].set_xticks(range(0,ncols,4), range(0,ncols,4))
                    ax[i].set_yticks(range(0,nrows,4), range(0,nrows,4));
                for i in [0,1,3,4]:
                    ax[i].legend()

                fig.suptitle(reference_gene + ' vs ' + genes[gidx] + ' : ' + bname + ' : Difference distribution', fontsize=fs+3)
                fig.tight_layout()

                filename = ddst + 'difference_detection_summary_' + bname + '.jpg'
                print('Saved image',filename)
                plt.savefig(filename, format='jpg', dpi=72, bbox_inches='tight', pil_kwargs={'optimize':True})
                plt.close()

if __name__ == '__main__':
    main()
