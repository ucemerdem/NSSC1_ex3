import matplotlib.pyplot as plt

threads = [1, 5, 10, 20, 40, 80]
sinx= [0.81533, 0.282098, 0.198622, 0.127703, 0.0704629, 0.0600694]
cos2xinv = [1.01566, 0.359375, 0.242065, 0.144327, 0.081681, 0.0587195]
x4m5= [0.596725, 0.251936, 0.164719, 0.109085, 0.0539882, 0.0437341]

def calc_speedup(func):
    listi = []
    thread_1 = func[0]
    thread_p = func[1:]
    for item in thread_p:
        listi.append(thread_1/item)
    return listi


def plot_runtime_threads(threads, sinx, cos2xinv, x4m5):
    plt.figure(figsize=(8, 6))
    plt.plot(threads, x4m5, label="X4M5", marker="o", color="deepskyblue")
    plt.plot(threads, sinx, label="SINX", marker="o", color="firebrick")
    plt.plot(threads, cos2xinv, label="COS2XINV", marker="o", color="darkorange")
    plt.xlabel("Threads Used")
    plt.ylabel("Runtime (s)")
    plt.title("Fun with the amazing HPC Cluster")
    plt.legend()
    plt.tight_layout()
    plt.show()

def plot_speed_up(threads, sinx, cos2xinv, x4m5):
    speedup_sin = calc_speedup(sinx)
    speedup_cos = calc_speedup(cos2xinv)
    speedup_x4 = calc_speedup(x4m5)
    plt.figure(figsize=(8, 6))
    plt.plot(threads[1:], speedup_x4, label="X4M5", marker="o", color="deepskyblue")
    plt.plot(threads[1:], speedup_sin, label="SINX", marker="o", color="firebrick")
    plt.plot(threads[1:], speedup_cos, label="COS2XINV", marker="o", color="darkorange")
    plt.xlabel("Threads Used")
    plt.ylabel("Speedup")
    plt.title("A beutiful Speedup Plot")
    plt.legend()
    plt.tight_layout()
    plt.show()

plot_runtime_threads(threads, sinx, cos2xinv, x4m5)
plot_speed_up(threads, sinx, cos2xinv, x4m5)

