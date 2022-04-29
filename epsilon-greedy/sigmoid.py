import matplotlib.pyplot as plt
import math

def sigmoid(x, h):
    return -0.5+1/(1+math.exp(-h*x))

xs = [i/10 for i in range(-100, 100)]
ys = [sigmoid(xs[i], 1) for i in range(len(xs))]

print(sigmoid(0, 1))

plt.plot(xs, ys)
plt.show()