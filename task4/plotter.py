import matplotlib.pyplot as plt

threads = [1, 2, 4, 8, 10, 20, 40]
static = [1607.12, 1383.61, 1410.47, 796.474, 574.231, 223.734, 173.804]
static1 = [1062.05, 1074.59, 1338.26, 1116.32, 815.824, 359.266, 249.827]
dynamic= [758.044, 676.84, 666.89, 523.176, 510.006, 277.205, 225.815]

def calc_speedup(func):
    listi = []
    thread_1 = func[0]
    thread_p = func[1:]
    for item in thread_p:
        listi.append(thread_1/item)
    return listi


def plot_runtime_threads(threads, static, static1, dynamic):
    plt.figure(figsize=(8, 6))
    plt.plot(threads, static, label="static", marker="o", color="deepskyblue")
    plt.plot(threads, static1, label="static,1", marker="o", color="firebrick")
    plt.plot(threads, dynamic, label="dynamic", marker="o", color="darkorange")
    plt.xlabel("Threads Used")
    plt.ylabel("Runtime (s)")
    plt.title("A not so wild Runtime Plot anymore")
    plt.legend()
    plt.tight_layout()
    plt.show()

def plot_speed_up(threads, static, static1, dynamic):
    speedup_s = calc_speedup(static)
    speedup_s1 = calc_speedup(static1)
    speedup_d = calc_speedup(dynamic)
    plt.figure(figsize=(8, 6))
    plt.plot(threads[1:], speedup_s, label="static", marker="o", color="deepskyblue")
    plt.plot(threads[1:], speedup_s1, label="static,1", marker="o", color="firebrick")
    plt.plot(threads[1:], speedup_d, label="dynamic", marker="o", color="darkorange")
    plt.xlabel("Threads Used")
    plt.ylabel("Speedup")
    plt.title("A beautiful Speedup Plot")
    plt.legend()
    plt.tight_layout()
    plt.show()

plot_runtime_threads(threads, static, static1, dynamic)
plot_speed_up(threads, static, static1, dynamic)

