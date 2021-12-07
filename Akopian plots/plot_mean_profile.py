import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import os
import sys
# print "This is the name of the script: ", sys.argv[0]
import locale

locale.setlocale(locale.LC_ALL, '')


# print(locale.format("%e", 3.14))
# print(locale.format("%f", 3.14))

def main(argv):
    try:
        directory = argv[0]
    except:
        directory = "output/"
    if directory[-1] != '/':
        directory += '/'
    if not os.path.exists(directory):
        print('ERROR: Output directory not found.')
        raise

    try:
        saveto = argv[1]
    except:
        saveto = "images_mean/"
    if saveto[-1] != '/':
        saveto += '/'
    if not os.path.exists(saveto):
        os.makedirs(saveto)

    files = os.listdir(directory)
    files.sort()

    # aaa = 0
    # print(aaa123456)

    f = open('text.txt', 'w')

    ccc = 0
    for file in files:
        if file.endswith(".dat") and '_' in file:
            ind = file.index('phi_')
            phi = float(file[ind + 4:len(file) - 4])
            dat = np.loadtxt(directory + file, delimiter=',')
            sz = int(np.sqrt(dat.shape[0]))
            dat = dat.reshape(sz, sz, 3)

            # L = [1, 2, 3, 4, 5]
            # aaa = 0

            aaa = 0
            bbb = 0
            for i in range(len(dat)):
                for j in range(len(dat)):
                    aaa = aaa + 1
                    if dat[i, j, 2] == dat[i, j, 2]:
                        bbb = bbb + dat[i, j, 2]
            #               print(dat[i], end=' ')

            #                aaa = aaa + 1

            #            f.write(ccc, bbb, '\n')
            #            f.write(str(ccc) + ' ' + str(bbb) + '\n')
            #            f.write(str(int(bbb)) + '\n')

            # print(locale.format("%e", 3.14))

            f.write('{0}\n'.format(str(locale.format("%e", bbb))))
            # f.write(str(ccc) + '\t'  + str(locale.format("%e", bbb)) + '\n')
            # f.write(str(bbb) + '\n')
            ccc = ccc + 1
            print("123", aaa, "     ", bbb)
            # print("123", dat[0,1,2])
            fig = plt.figure(figsize=(6, 6))
            ax = plt.subplot(111)
            ax.imshow(dat[:, :, 2], origin='lower',
                      extent=(dat[:, :, 0].min(), dat[:, :, 0].max(), dat[:, :, 1].min(), dat[:, :, 1].max()))
            ax.set_title(r'$\phi={{{}}}^{{\circ}}$'.format(phi))
            ax.set_xlabel(r'$a_2 / R_{\rm star}$')
            ax.set_ylabel(r'$a_1 / R_{\rm star}$')
            plt.savefig(saveto + file[:-4] + '.png')
            plt.close()
    f.close()


if __name__ == '__main__':
    main(sys.argv[1:])
