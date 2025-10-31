import sys
import pandas as pd
import matplotlib.pyplot as plt

sep = '\t'

args = sys.argv

data_path = f'../results/{args[1]}.txt'
plot_path = f'../plots/{args[1]}.png'

y_max = 30

with open(data_path) as header:
    title, x_label, y_label = header.readline().split(sep)

data = pd.read_csv(data_path, sep=sep, skiprows=1)

for i in range(1, data.shape[1]):
    plt.plot(data.iloc[:, 0], data.iloc[:, i])

plt.ylim(0, y_max)

plt.grid()
plt.xlabel(x_label)
plt.ylabel(y_label)
plt.title(title)

plt.savefig(plot_path, dpi=300)
