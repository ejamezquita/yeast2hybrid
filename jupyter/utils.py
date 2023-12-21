from PIL import Image
import numpy as np

# I assume that the picture has a green tape on of its sides
# Change `color_check = 2` if tape is blue

def load_image(filename, color_check=1, check_rotation=False):
    with open(filename, "rb") as fp:
        im = Image.open(fp)
        im.load()
    rgb = np.asarray(im.convert('RGB'))
    
    if not check_rotation:
        return rgb
        
    border = 300

    if im.height > im.width:
        topg = rgb[:border, :, color_check]
        botg = rgb[-border:, :, color_check]
        rot = 90
    else:
        topg = rgb[:,:border, color_check]
        botg = rgb[:, -border:, color_check]
        rot = 0
    if np.quantile(topg, .9) > np.quantile(botg, .9) :
        im = im.rotate(rot, expand=True)
    else:
        im = im.rotate(rot+180, expand=True)
        
    rgb = np.asarray(im.convert('RGB'))
    
    return rgb
