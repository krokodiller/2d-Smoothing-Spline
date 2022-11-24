import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
mpl.use('TkAgg')  # or can use 'TkAgg', whatever you have/prefer


def print_heatmap():
    v = []
    file = open("solution.txt")
    data = file.readline().split(" ")
    left = float(data[0])
    bottom = float(data[1])
    right = float(data[2])
    top = float(data[3])
    for line in file:
        data = line.split()
        vline = []
        for value in data:
            vline.append(float(value))
        v.append(vline)
    file.close()

    vsidecount = len(v)
    xl = np.linspace(left, right, vsidecount)
    yl = np.linspace(bottom, top, vsidecount)
    x, y = np.meshgrid(xl, yl)

    plt.contourf(x, y, v, 100)
    plt.show()

def print_surface():
    file = open("solution.txt")


    data = file.readline().split(" ")
    left = float(data[0])
    bottom = float(data[1])
    right = float(data[2])
    top = float(data[3])
    file.close()

    v = np.genfromtxt('solution.txt', skip_header=1)

    
    size = np.shape(v)[0]
    xl = np.linspace(left, right, size)
    yl = np.linspace(bottom, top, size)
    x, y = np.meshgrid(xl, yl)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    surf = ax.plot_surface(x, y, v, cmap=cm.coolwarm, linewidth=0)

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('F')
    fig.colorbar(surf, shrink=0.5, aspect=5)

    plt.show()



if __name__ == "__main__":
    print_surface()