import math
from matplotlib import pyplot as plt

def healthfn(x):
    if x <= 20:
        return 2*math.log(10*x)
    else:
        return 2*(-(53/64000)*(x**2) + (53/1600)*x + (795/160))

xs = [i/10 for i in range(1, 1001)]
ys = [healthfn(x) for x in xs]

plt.plot(xs, ys)
plt.show()

# quadratic with f(20) = 52.98317 and f(100) = 0