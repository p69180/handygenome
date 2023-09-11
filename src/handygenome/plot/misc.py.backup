import colorsys

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt



def lightness_spectrum(color, n, l_start=0.3, l_end=0.7):
    color_hls = colorsys.rgb_to_hls(*mpl.colors.to_rgb(color))
    return [
        colorsys.hls_to_rgb(color_hls[0], x, color_hls[2]) 
        for x in np.linspace(l_start, l_end, n + 2, endpoint=True)[1:-1]
    ]

    
