import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 20})
plt.figure(figsize=(12, 8))
args = sys.argv

data_path = f'../results/{args[1]}.csv'
plot_path = f'../plots/{args[1]}.png'

with open(data_path) as header:
    title, x_label, y_label = header.readline().rstrip().split('\t')
    y_min, y_max = map(float, header.readline().rstrip().split('\t'))
    legend_labels = header.readline().rstrip().split('\t')

has_legend = not (len(legend_labels) == 1 and legend_labels[0] == '')

data = pd.read_csv(data_path, skiprows=3, header=None)

ls = ['-', '-', '--', '--']
cs = ['C0', 'C1', 'C0', 'C1']

for i in range(1, data.shape[1]):
    # plt.plot(data.iloc[:, 0], data.iloc[:, i], label='' if not has_legend else legend_labels[i - 1], ls=ls[i - 1], lw=5, zorder=-i, c=cs[i - 1])
    plt.scatter(data.iloc[:, 0], data.iloc[:, i], s=3, label='' if not has_legend else legend_labels[i - 1], zorder=-i)

# plt.axhline(0.285, ls='--', c='black', lw=3)
# plt.text(
#     x=-1.3,
#     y=0.285,
#     s=rf'$1.5 \lambda$',
#     color='black',
#     va='center',
#     ha='right'
# )

y_lim = y_min != 0 or y_max != 0
y_min = y_min if y_lim else None
y_max = y_max if y_lim else None
plt.ylim(y_min, y_max)

# plt.ticklabel_format(axis='x', useOffset=data.iloc[0, 0])

plt.minorticks_on()
plt.grid(which='major', color='#666666', linestyle='-', alpha=0.6)
plt.grid(which='minor', color='#999999', linestyle=':', alpha=0.6)

plt.xlabel(x_label)
plt.ylabel(y_label)
# plt.title(title)

# if has_legend: 
#     plt.legend(loc='upper left') 
plt.savefig(plot_path, dpi=300)
