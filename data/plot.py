import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

data1 = pd.read_csv('/Users/xzoidy/MFF/Nuclear_Computing/data/6Li_Dytrich.txt', sep='\s+', header=None)
data1 = pd.DataFrame(data1)

Nmax = data1[0]
hw = data1[1]
E = data1[2]

plt.xlabel(r'$h\omega$')
plt.ylabel(r'E')

# Use a dictionary to map Nmax values to unique colors
nmax_color_map = {2: "blue", 4: "green", 6: "red"}

# Group data by Nmax
grouped_data = data1.groupby(0)

for nmax_value, group in grouped_data:
    color = nmax_color_map[nmax_value]
    
    # Plot points
    plt.scatter(group[1], group[2], color=color, label=f'Nmax = {nmax_value}')
    
    # Connect points with lines
    plt.plot(group[1], group[2], color=color, linestyle='-', marker='o', label='')

plt.legend()
plt.show()