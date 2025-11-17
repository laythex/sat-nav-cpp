import sys
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

args = sys.argv

data_path = f'../results/{args[1]}.csv'

with open(data_path) as header:
    title, mode = header.readline().rstrip().split('\t')
    threshold = float(header.readline().rstrip())

plot_path = f'../plots/{args[1]}-{mode}.png'

data = pd.read_csv(data_path, skiprows=2)

plt.figure(figsize=(7.5, 4))

m = Basemap(projection='cyl', resolution='c')
m.fillcontinents(color='gray',lake_color='white')
m.drawparallels(range(-90, 91, 30))
m.drawmeridians(range(-180, 181, 60))

if mode == 'full':
    plt.scatter(data.iloc[:, 0], data.iloc[:, 1], c=data.iloc[:, 2], s=3)
    plt.colorbar()
elif mode == 'threshold':
    mask = data.iloc[:, 2] > threshold
    plt.scatter(data[mask].iloc[:, 0], data[mask].iloc[:, 1], s=3)

plt.title(f'{title} > { threshold} м' if mode == 'threshold' else title)

plt.savefig(plot_path, dpi=300)