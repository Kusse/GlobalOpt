import numpy as np
import matplotlib.pyplot as plt

imfile = "sim_im55.txt"
im = np.loadtxt(imfile)
plt.imshow(im)
#plt.xlim([-20,20])
#plt.ylim([-20,20])
#plt.pcolor(im) 
#plt.axis("equal")
plt.show() 
    