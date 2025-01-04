import matplotlib.pyplot as plt

threads = [1, 2, 4, 8, 10, 20, 40]
static = [636.512, 1076.68, 1117.12, 631.589, 459.566, 433.039, 311.395]
static1 = [657.518, 1005.71, 1216.18, 989.14, 687.595, 603.176, 398.109]
dynamic= [801.76, 686.001, 676.06, 448.696, 472.041, 380.726, 284.526]

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
    plt.title("A very wild Runtime Plot D:")
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
    plt.title("A strange Speedup Plot")
    plt.legend()
    plt.tight_layout()
    plt.show()

plot_runtime_threads(threads, static, static1, dynamic)
plot_speed_up(threads, static, static1, dynamic)

