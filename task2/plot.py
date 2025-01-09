import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import subprocess

def plot_and_save_residuals(filename, x, y_residuals):
    # this basically holds the number of iterations carried out for each approach of the CG algorithm
    xrange = len(x)

    fig, ax = plt.subplots()
    ax.set(xlabel="CG iterations", ylabel=r"$\frac{{\| r_k \|}_2}{{\| r_0 \|}_2}$",title=f'{filename}')

    plt.plot(x, y_residuals[:xrange], '-', color="red", label="Own Implementation")
    plt.plot(x, y_residuals[xrange:2*xrange], '-', color="blue", label="Eigen DiagonalPreconditioner")
    plt.plot(x, y_residuals[2*xrange:], '-', color="orange", label="Eigen IncompleteCholesky")
    plt.yscale("log")
    plt.legend()

    fig.savefig(filename + "_residuals_" + str(len(x)) + ".png")

    plt.show()


def run_c_program():
    filename = "bcsstk11.mtx"
    iters = 100
    norms_residuals = []

    args = [filename, str(iters)]
    # Run the C++ program and capture stdout
    result = subprocess.run(
        [os.path.join(sys.path[0],"cg")] + args, 
        capture_output=True,     # Capture stdout
        text=True                # Decode stdout to string
    )
    # C++ program output
    output = result.stdout

    for line in output.splitlines():
        norms_residuals.append(np.float64(line))
        # this list contains all residuals from the three different approaches of the CG algorithm (own implementation, Eigen DiagonalPreconditioner, Eigen IncompleteCholesky)
        # the first 'iters' elements are the residuals from the own implementation
        # the second 'iters' elements are the residuals from the Eigen DiagonalPreconditioner
        # the third 'iters' elements are the residuals from the Eigen Incomplete

    plot_and_save_residuals(args[0]+"_performance", np.arange(iters), norms_residuals)


run_c_program()
