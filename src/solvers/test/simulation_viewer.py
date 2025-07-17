import sympy as sym
import numpy as np

# Plotting stuff
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
from matplotlib.ticker import MaxNLocator
import copy
#matplotlib.use('QT5Agg') # Or 'QT4Agg'

params = {'legend.fontsize': 18,
        'axes.labelsize': 18,
        'axes.titlesize': 18,
        'xtick.labelsize': 18,
        'ytick.labelsize': 18}
plt.rcParams.update(params)
#plt.rcParams['figure.figsize'] = [5, 2.5] # for square canvas

def main():

    #main_1D()
    main_2D()

def main_2D():

    show_man = not False

    nTimeSteps = 7
    n = 24   # n DOFs for each variable
    N = n**2 # Total # of DOFs
    data_dir = "../../../build/debug/src/solvers/test/out/n" + str(n) + "/"

    ax_lims_x = []
    ax_lims_y = []

    ax_lims_x2= []
    ax_lims_y2 = []

    x = np.linspace(0, 1, n)
    y = x
    X, Y = np.meshgrid(x, y)

    for step in range(0, nTimeSteps+1, 1):

        uman_file = data_dir + "uman_" + str(step) + ".ascii"
        unum_file = data_dir + "unum_" + str(step) + ".ascii"        

        with open(uman_file) as f:
            header = f.readline()
            import re
            time = re.findall(r'"(.*?)"', header)[0]
        
        uman_data = np.genfromtxt(uman_file, skip_header=1)
        unum_data = np.genfromtxt(unum_file, skip_header=1)

        uman_data = uman_data[0:N].reshape(n, n)
        unum_data = unum_data[0:N].reshape(n, n)

        # Manufactured
        fig2, ax2  = plt.subplots(subplot_kw={"projection": "3d"})
        fig2.canvas.manager.window.move(200, 200)
        col_map = copy.copy(cm.coolwarm)
        surf = ax2.plot_surface(X, Y, uman_data - unum_data, cmap=col_map)
        ax2.view_init(azim = 225)
        ax2.set_xlabel("$x$")
        ax2.set_ylabel("$y$")
        ax2.grid()
        ax2.set_title( "manufactured: time = {}, step = {}".format( time, step ) )
        if step == 0:
            ax_lims_x2 = ax2.get_xlim()
            ax_lims_y2 = ax2.get_ylim()
        else:
            ax2.set_xlim( ax_lims_x2 )
            ax2.set_ylim( ax_lims_y2 )

        # Numerical
        fig, ax  = plt.subplots(subplot_kw={"projection": "3d"})
        fig.canvas.manager.window.move(1000, 200)
        col_map = copy.copy(cm.coolwarm)
        surf = ax.plot_surface(X, Y, unum_data, cmap=col_map)
        ax.view_init(azim = 225)
        ax.set_xlabel("$x$")
        ax.set_ylabel("$y$")
        ax.grid()
        ax.set_title( "numerical: time = {}, step = {}".format( time, step ) )
        if step == 0:
            ax_lims_x = ax.get_xlim()
            ax_lims_y = ax.get_ylim()
        else:
            ax.set_xlim( ax_lims_x )
            ax.set_ylim( ax_lims_y )

        

        plt.show()

def main_1D():

    show_man = not False

    nTimeSteps = 300
    n = 128 # n-1 DOFs for each variable
    data_dir = "../../../build/debug/src/solvers/test/out/n" + str(n) + "/"

    ax_lims_x = []
    ax_lims_y = []

    for step in range(0, nTimeSteps+1, 10):

        uman_file = data_dir + "uman_" + str(step) + ".ascii"
        unum_file = data_dir + "unum_" + str(step) + ".ascii"        

        with open(uman_file) as f:
            header = f.readline()
            import re
            time = re.findall(r'"(.*?)"', header)[0]
        
        uman_data = np.genfromtxt(uman_file, skip_header=1)
        unum_data = np.genfromtxt(unum_file, skip_header=1)

        uman_data = uman_data[0:n-1]
        unum_data = unum_data[0:n-1]

        fig, ax = plt.subplots()

    # Example: Move the window to a specific position (e.g., 100 pixels from top, 200 pixels from left)
        fig.canvas.manager.window.move(200, 200)

        if show_man:
            ax.plot(uman_data, '--bo', markersize = 10, markerfacecolor = "white", label = "$u_{{man}}$")
        ax.plot(unum_data, '-r>', label = "$u_{{num}}$")

        #ax.plot(unum_data - uman_data, '-r>', label = "$e_{{num}}$")
            

        # print(Tman_data, sep = "\n\n")
        # print(Tnum_data)

        ax.set_xlabel("$x$")
        ax.legend()
        ax.grid()
        ax.set_title( "time = {}, step = {}".format( time, step ) )

        if step == 0:
            ax_lims_x = ax.get_xlim()
            ax_lims_y = ax.get_ylim()
        else:
            ax.set_xlim( ax_lims_x )
            ax.set_ylim( ax_lims_y )

        plt.show()
 


# Call the main method!
if __name__ == "__main__":
    main()