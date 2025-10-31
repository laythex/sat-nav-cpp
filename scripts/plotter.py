import sys
import pandas as pd
import matplotlib.pyplot as plt

args = sys.argv

data_path = f'../results/{args[1]}.csv'
plot_path = f'../plots/{args[1]}.png'

with open(data_path) as header:
    title, x_label, y_label = header.readline().split('\t')
    y_min, y_max = map(float, header.readline().split('\t'))

data = pd.read_csv(data_path, skiprows=2)

for i in range(1, data.shape[1]):
    plt.plot(data.iloc[:, 0], data.iloc[:, i])

plt.ylim(y_min, y_max)

plt.grid()
plt.xlabel(x_label)
plt.ylabel(y_label)
plt.title(title)

plt.savefig(plot_path, dpi=300)
