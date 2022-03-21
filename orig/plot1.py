import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import os
import sys
# print "This is the name of the script: ", sys.argv[0]

def main(argv):
    try:
        directory = argv[0]
    except:
        directory = "output/"
    if (directory[-1] != '/'):
        directory += '/'
    if not os.path.exists(directory):
        print ('ERROR: Output directory not found.')
        raise

    try:
        saveto = argv[1]
    except:
        saveto = "images/"
    if (saveto[-1] != '/'):
        saveto += '/'
    if not os.path.exists(saveto):
        os.makedirs(saveto)

    files = os.listdir(directory)
    files.sort()

    for file in files:
        if file.endswith(".dat") and '_' in file:
            ind = file.index('phi_')
            phi = float(file[ind+4:len(file)-4])
            dat = np.loadtxt(directory + file, delimiter=',')
            sz = int(np.sqrt(dat.shape[0]))
            dat = dat.reshape(sz, sz, 3)

            fig = plt.figure(figsize=(6, 6))
            ax = plt.subplot(111)
            ax.imshow(dat[:,:,2], origin='lower', extent=(dat[:,:,0].min(), dat[:,:,0].max(), dat[:,:,1].min(), dat[:,:,1].max()))
            ax.set_title(r'$\phi={{{}}}^{{\circ}}$'.format(phi))
            ax.set_xlabel(r'$a_2 / R_{\rm star}$')
            ax.set_ylabel(r'$a_1 / R_{\rm star}$')
            plt.savefig(saveto + file[:-4] + '.png')
            plt.close();

if __name__ == '__main__':
    main(sys.argv[1:])
