import sys
import pandas as pd
import matplotlib.pyplot as plt

args = sys.argv

data_path = '../results/' + args[1]
plot_path = '../plots/' + args[2]

y_max = 75

with open(data_path) as header:
    title = header.readline()

data = pd.read_csv(data_path, sep='\t', skiprows=1)
x_label, y_label = data.columns.tolist()

plt.plot(data[x_label], data[y_label])

plt.ylim(0, y_max)

plt.grid()
plt.xlabel(x_label)
plt.ylabel(y_label)
plt.title(title)

plt.savefig(plot_path, dpi=300)
