# Yeast2Hybrid: Image processing to identify color differences

**Author**: Erik Amézquita

**Date**: December 2023

**Please read carefully before start**

## Description

This is a collection of 4 jupyter notebooks to process yeast plate images, extract color intensities for each yeast colony, compare those color intensities with reference plates, and identify those colonies with the highest intensity differences.

Disclaimer: I am **not** a software developer. This is **not** meant to be a plug-and-play package for all kind of yeast image processing. It is rather a breakdown of all the steps needed to process a set of particular images. Please read carefully each of the jupyter notebooks and adapt them to your data. I tried my best to explain all my steps and point out which bits could potentially change with a different set of images.

The first part (the first jupyter notebook) is meant to be manually run for each yeast image. You are encourage to toggle several image processing parameters that eventually yield the centers of each yeast colony. All these manually set parameters are then automatically saved in a metadata file. 

In principle, with this metadata, the following parts can be run automatically as a `py` script. You are however encouraged to read the corresponding notebooks to make sure they work as intended with *your* own data. Feel free to adapt the notebooks to your needs and then remember to change the `py` files accordingly.

## Python packages

All the code here uses common python libraries that are easily installable.

- Numpy
- Matplotlib
- Scipy
- Pandas
- Pillow
- Skimage

All these libraries come preinstalled with anaconda. The rest of used libraries (os, argparse, glob, importlib) should come with any python distribution.

## Work directory setup

The notebooks assume that your work directory follows a specific structure. Folder and filenames can change as long as they follow the specific structure as in this repository. 

- Make sure all your files are in a single directory (I'll refer to this directory as the *root* directory).
- Inside this root, **make sure to have 4 directories**.
    - `diagnostic`: where summary images will be created for final visual assessment of the correctness of the pipeline
    - `jupyter`: where all the jupyter notebooks and py files are
    - `proc`: processed data and numerical summaries
    - `raw`: where all raw images to be analyzed are located

- In root you should have a CSV with the reference gene information.
    
- In principle, `diagnostic` and `proc` are empty folders at first, where adecuate subfolders and files are automatically stored.

### Structure of the `raw` folder

- Inside `raw`, there should be separate folders for different sources of images.
    - In this case, you'll see `mather` and `leyre` subfolders refering to images taken by different colleagues under different conditions.
    - `mather` files are our reference or control
    - `leyre` files represent mutations
    
- Inside each of source subfolder, you'll have one subsubfolder per gene to analyze
    - In this case, `leyre`'s images are related to gene *LBD37*
    
- Inside each gene subsubfolder, you'll have the actual raw images
    - The actual name of the images is not important as long as they are consistent
    - **Make sure** that somewhere in the filename they say `_plate_XX_` to indicate which plate number we are looking at
    - Make sure that all folders have the same number of plate images
    - Additional raw images (like a color checker) are ok
    - **Do not** place two `plate_01` images in the same subsubfolder
    
- You'll see that this structure is then repeated for the `proc` and `diagnostic` folders


### The jupyter notebooks

- The pipeline is split into 4 parts.
- It is assumed that all your data has undergone Part I, before it can processed in Part II
- Same idea for Parts III and IV
- The jupyter notebook filenames should make obvious which part is which
- I tried my best to write plenty of useful comments in the notebooks
- Please read carefully the notebooks for further details

#### Part I: Find colony centers

- **Input** 
    - A raw image of yeast plate
- **Output**
    - A CSV with the image processing parameters used
    - A NPY file (numpy's binary save file) with the (X,Y)-coordinates of the colony centers
    
![](https://raw.githubusercontent.com/ejamezquita/ejamezquita.github.io/main/yeast2hybrid/figs/raw_plate_5.jpg)

![](https://raw.githubusercontent.com/ejamezquita/ejamezquita.github.io/main/yeast2hybrid/figs/centers_plate_5.jpg)

- This notebook is meant to be run manually for every single image.
- If all your raw images were taken under the same conditions, then the needed image processing parameters should be the same or quite similar.
- Choosing the right image processing parameters will ensure the correctness of the ensuing pipeline

#### Part II: Compute color intensity

- **Input**
   - A raw image of yeast plate
   - A CSV with the image processing parameters used
   - A NPY file (numpy's binary save file) with the (X,Y)-coordinates of the colony centers
- **Output**
   - A CSV with the median color intensity for each colony after correcting for nonuniform illumation
   - (Optional) A diagnostic image which shows illumination correction, background removal, and the final intensity matrix
   
![](https://raw.githubusercontent.com/ejamezquita/ejamezquita.github.io/main/yeast2hybrid/figs/light_correction_20231130_plate_05.jpg)

- This part ideally is meant to be run as a `py` script
- The script expects parameters that specify the workdirectory structure
- See the `02` py file for more details
- The py script first checks whether the output CSV file already exists
    - If it doesn't exist, it computes the CSV
    - Else, it skips it and goes to the next image
- This way you don't have to recompute/overwrite the same data over and over again.

#### Part III: Compare color intensities

- **Input**
    - The intensity matrix from a mutant plate
    - The intensity matrix from the control/reference plate
    - A CSV with gene reference (what coordinate in the plate is what gene)
- **Output**
    - Two CSVs with the genes that exhibit the largest color intensity differences
    - There are two CSVs corresponding to two different criteria to define *large* intensity value difference.
    - (Optional) A diagnostic summary image which shows how intensities are distributed for each image, how these were matched, the distribution of intensity differences, and the criteria to determine what is *large* intensity difference.
    
![](https://github.com/ejamezquita/ejamezquita.github.io/blob/main/yeast2hybrid/figs/difference_detection_summary_20231130_plate_05.jpg)

- Same comments as in Part II

#### Part IV: Visualize candidates

- **Input**
   - A CSV with candidate genes of interest and their respective locations in the plate
   - The raw image of the mutant plate
   - The raw image of the reference plate
- **Output**
   - A diagnostic image that highlights the colonies that showcase the largest differences in color intensity.
   - These colonies are shown side by side in control and mutant plates for final visual assessment.
    
![](https://raw.githubusercontent.com/ejamezquita/ejamezquita.github.io/main/yeast2hybrid/figs/difference_detection_summary_20231130_plate_05.jpg)

- Same comments as in Part II
- Notice that colonies are arranged in descending order, with the one with the most reported color intensity difference first.

## What next?

- You might want to curate the CSVs from Part III by looking at Part IV
- You can also combine easily all CSVs in R or python for further wrangling and filtering.

