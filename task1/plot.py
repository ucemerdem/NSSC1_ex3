import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import subprocess

def plot_and_save_residuals(filename, x, y_residuals):
    fig, ax = plt.subplots()
    ax.set(xlabel="CG iterations", ylabel=r"$\frac{{\| r_k \|}_2}{{\| r_0 \|}_2}$",title=f'{filename}')
    plt.plot(x, y_residuals, '-', color="blue")
    fig.savefig(filename + "_residuals_" + str(len(y_residuals)) + ".png")
    plt.show()

def plot_and_save_errorsANorm(filename, x, y_A_norm):    
    fig, ax = plt.subplots()
    ax.set(xlabel="CG iterations", ylabel=r"error on A-norm $\| e_k \|$",title=f'{filename}')
    plt.plot(x, y_A_norm, '-', color="orange")
    fig.savefig(filename + "_errorA-norm_" + str(len(y_A_norm)) + ".png")
    plt.show()


def run_c_program():
    filename = "bcsstk13.mtx"
    iters = 100
    norms_residuals = []
    errors_A_norm = []

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
        norms_residuals.append(np.float64(str.split(line)[0]))
        errors_A_norm.append(np.float64(str.split(line)[1]))

    plot_and_save_residuals(args[0]+"_performance", np.arange(iters), norms_residuals)
    plot_and_save_errorsANorm(args[0]+"_performance", np.arange(iters), errors_A_norm)


run_c_program()
