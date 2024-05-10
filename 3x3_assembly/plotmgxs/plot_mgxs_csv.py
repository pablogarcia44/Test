import sys
import os 
import openmc
import numpy as np
import matplotlib.pyplot as plt
import openmc.mgxs as mgxs
import pandas as pd

#Choose reaction and cell to plot
reaction='fission'
cell='fuel'
s
csv_file_path_xs = os.path.join(os.path.dirname(os.path.abspath(__file__)), "mgxs", f"{reaction}_{cell}_mgxs.csv")
csv_file_path_energy = os.path.join(os.path.dirname(os.path.abspath(__file__)), "mgxs", f"energy_groups.csv")

df_xs = pd.read_csv(csv_file_path_xs)
df_energy = pd.read_csv(csv_file_path_energy)
mgxs = df_xs['mean'].values
energy = df_energy['0'].values
name_energy = df_energy['1'].values[1]


fig, ax = plt.subplots()
ax.step(energy[:-1], mgxs, where='post', label=name_energy) 
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('Energy [eV]')
ax.set_ylabel('MGXS [cm-1]')
#ax.set_xlim(10e-12, 10e-6)  
#ax.set_ylim(0.5*10e-3, 10e1)
ax.grid()
ax.legend()
plt.show()


