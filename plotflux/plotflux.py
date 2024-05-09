import sys
import os 
import openmc
import numpy as np
import matplotlib.pyplot as plt
import openmc.mgxs as mgxs
import math
import pickle
from IPython.display import Image
import openmc.mgxs


repertoire_courant = os.path.dirname(__file__)
chemin_xs = os.path.join(repertoire_courant, 'cross_sections.xml')

model=openmc.Model()

#Materials

uo2 = openmc.Material(name='uo2')
uo2.add_nuclide('U235', 0.03)
uo2.add_nuclide('U238', 0.97)
uo2.add_nuclide('O16', 2.0)
uo2.set_density('g/cm3', 10.0)

zirconium = openmc.Material(name="zirconium")
zirconium.add_element('Zr', 1.0)
zirconium.set_density('g/cm3', 6.6)

helium=openmc.Material(name='helium')
helium.add_element('helium',1)
helium.set_density('g/cm3',8e-4)

water = openmc.Material(name="h2o")
water.add_nuclide('H1', 2.0)
water.add_nuclide('O16', 1.0)
water.set_density('g/cm3', 1.0)
water.add_s_alpha_beta('c_H_in_H2O')

model.materials = openmc.Materials([uo2, zirconium, water,helium])

model.materials.cross_sections = chemin_xs

model.materials.export_to_xml()

#Geometry

#create surface
fuel_outer_radius = openmc.ZCylinder(r=0.39)
clad_inner_radius = openmc.ZCylinder(r=0.40)
clad_outer_radius = openmc.ZCylinder(r=0.46)

pitch = 1.26
left = openmc.XPlane(-pitch/2, boundary_type='reflective')
right = openmc.XPlane(pitch/2, boundary_type='reflective')
bottom = openmc.YPlane(-pitch/2, boundary_type='reflective')
top = openmc.YPlane(pitch/2, boundary_type='reflective')

#create region of interest
fuel_region = -fuel_outer_radius
gap_region = +fuel_outer_radius & -clad_inner_radius
clad_region = +clad_inner_radius & -clad_outer_radius
water_region = +left & -right & +bottom & -top & +clad_outer_radius

#create cell
fuel = openmc.Cell(name='fuel')
fuel.fill = uo2
fuel.region = fuel_region

gap = openmc.Cell(name='gap')
gap.fil=helium
gap.region = gap_region

clad = openmc.Cell(name='clad')
clad.fill = zirconium
clad.region = clad_region

moderator = openmc.Cell(name='moderator')
moderator.fill = water
moderator.region = water_region

root_universe = openmc.Universe(cells=(fuel, gap, clad, moderator))
model.geometry = openmc.Geometry(root_universe)

model.geometry.export_to_xml()

#Parameter

point = openmc.stats.Point((0, 0, 0))
source = openmc.Source(space=point)

batches = 50
inactive = 10
particles = 3000

settings = openmc.Settings()
settings = openmc.Settings()
settings.batches = batches
settings.inactive = inactive
settings.particles = particles
##settings.output = {'tallies': False}
#settings.run_mode = 'eigenvalue'
#settings.verbosity = 4

point = openmc.stats.Point((0, 0, 0))
source = openmc.Source(space=point)

model.settings = settings

model.settings.export_to_xml()

#MGXS


#create energy group discretisation
VIMS69=np.array([1e-11, 5e-09, 1e-08, 1.5e-08, 2e-08, 2.5e-08, 3e-08, 3.5e-08, 4.2e-08, 5e-08, 5.8e-08, 6.7e-08, 8e-08, 1e-07, 1.4e-07, 1.8e-07, 2.2e-07, 2.5e-07, 2.8e-07, 3e-07, 3.2e-07, 3.5e-07, 4e-07, 5e-07, 6.25e-07, 7.8e-07, 8.5e-07, 9.1e-07, 9.5e-07, 9.72e-07, 9.96e-07, 1.02e-06, 1.045e-06, 1.071e-06, 1.097e-06, 1.123e-06, 1.15e-06, 1.3e-06, 1.5e-06, 2.1e-06, 2.6e-06, 3.3e-06, 4e-06, 9.877e-06, 1.5968e-05, 2.77e-05, 4.8052e-05, 7.55014e-05, 0.000148729, 0.000367263, 0.000906899, 0.0014251, 0.00223945, 0.0035191, 0.00553, 0.009118, 0.01503, 0.02478, 0.04085, 0.06734, 0.111, 0.183, 0.3025, 0.5, 0.821, 1.353, 2.231, 3.679, 6.0655, 10.0])

#openmc.mgxs.GROUP_STRUCTURES.keys() contain XMAS172
XMAS172=np.array([1.00001e-11, 3e-09, 5e-09, 6.9e-09, 1e-08, 1.5e-08, 2e-08, 2.5e-08, 3e-08, 3.5e-08, 4.2e-08, 5e-08, 5.8e-08, 6.7e-08, 7.7e-08, 8e-08, 9.5e-08, 1e-07, 1.15e-07, 1.34e-07, 1.4e-07, 1.6e-07, 1.8e-07, 1.89e-07, 2.2e-07, 2.48e-07, 2.8e-07, 3e-07, 3.145e-07, 3.2e-07, 3.5e-07, 3.91e-07, 4e-07, 4.33e-07, 4.85e-07, 5e-07, 5.4e-07, 6.25e-07, 7.05e-07, 7.8e-07, 7.9e-07, 8.5e-07, 8.6e-07, 9.1e-07, 9.3e-07, 9.5e-07, 9.72e-07, 9.86e-07, 9.96e-07, 1.02e-06, 1.035e-06, 1.045e-06, 1.071e-06, 1.097e-06, 1.11e-06, 1.123e-06, 1.15e-06, 1.17e-06, 1.235e-06, 1.3e-06, 1.3375e-06, 1.37e-06, 1.44e-06, 1.475e-06, 1.5e-06, 1.59e-06, 1.67e-06, 1.755e-06, 1.84e-06, 1.93e-06, 2.02e-06, 2.1e-06, 2.13e-06, 2.36e-06, 2.55e-06, 2.6e-06, 2.72e-06, 2.76792e-06, 3.3e-06, 3.38075e-06, 4e-06, 4.12925e-06, 5.043477e-06, 5.34643e-06, 6.160116e-06, 7.523983e-06, 8.315287e-06, 9.189814e-06, 9.905554e-06, 1.122446e-05, 1.370959e-05, 1.592827e-05, 1.945484e-05, 2.260329e-05, 2.49805e-05, 2.760773e-05, 3.051126e-05, 3.372015e-05, 3.726653e-05, 4.0169e-05, 4.551744e-05, 4.82516e-05, 5.157802e-05, 5.559513e-05, 6.790405e-05, 7.567357e-05, 9.166088e-05, 0.000136742, 0.0001486254, 0.000203995, 0.0003043248, 0.0003717032, 0.0004539993, 0.0006772874, 0.0007485183, 0.0009142423, 0.001010394, 0.001234098, 0.001433817, 0.001507331, 0.002034684, 0.002248673, 0.003354626, 0.003526622, 0.005004514, 0.005530844, 0.007465858, 0.00911882, 0.01113775, 0.01503439, 0.01661557, 0.02478752, 0.02739445, 0.029283, 0.03697864, 0.04086771, 0.05516564, 0.06737947, 0.08229747, 0.11109, 0.1227734, 0.1831564, 0.2472353, 0.2732372, 0.3019738, 0.407622, 0.450492, 0.4978707, 0.5502322, 0.6081006, 0.82085, 0.9071795, 1.002588, 1.108032, 1.224564, 1.353353, 1.652989, 2.018965, 2.231302, 2.46597, 3.011942, 3.678794, 4.49329, 5.488116, 6.065307, 6.7032, 8.187308, 10.0, 11.618343, 13.840307, 14.918247, 17.33253, 19.64033])



energies_xmas = openmc.mgxs.GROUP_STRUCTURES['XMAS-172']
xmas_filter = openmc.EnergyFilter(openmc.mgxs.GROUP_STRUCTURES['XMAS-172'])
tally = openmc.Tally()
tally.filters = [xmas_filter]
tally.scores = ['flux']
model.tallies = [tally]

statepoint_filename = model.run()
ce_spfile = './statepoint_ce.h5'
os.rename(statepoint_filename, ce_spfile)
ce_sumfile = './summary_ce.h5'
os.rename('summary.h5', ce_sumfile)
sp = openmc.StatePoint(ce_spfile, autolink=False)
su = openmc.Summary(ce_sumfile)
sp.link_with_summary(su)


t = sp.tallies[tally.id]
flux_mean = t.mean.ravel()
flux_unc = t.std_dev.ravel()
    
fig, ax = plt.subplots()
ax.step(energies_xmas[:-1], flux_mean/np.diff(energies_xmas), where='post', label='XMAS-172')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('Energy [eV]')
ax.set_ylabel('Flux [n-cm/eV-src]')
ax.grid()
ax.legend()
plt.show()





















 






