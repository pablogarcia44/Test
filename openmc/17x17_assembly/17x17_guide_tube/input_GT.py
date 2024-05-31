import sys
import os 
import openmc
import numpy as np
import matplotlib.pyplot as plt
import openmc.mgxs as mgxs
import pandas as pd
from IPython.display import Image
import glob

clean_file_list=["model.xml","materials.xml","geometry.xml","plots.xml","geometry.xml","settings.xml","tallies.out","summary.h5","statepoint.50.h5"]
for file in clean_file_list :
    path_file = os.path.join(os.getcwd(), file)
    if os.path.exists(path_file):
        os.remove(path_file)

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
water = openmc.Material(name="water")
water.add_nuclide('H1', 2.0)
water.add_nuclide('O16', 1.0)
water.set_density('g/cm3', 1.0)
water.add_s_alpha_beta('c_H_in_H2O')
model.materials = openmc.Materials([uo2, zirconium, water]) 
#Perso path
#model.materials.cross_sections = '/home/pablo/internship/xs_for_input/cross_sections.xml'
#Mac path
model.materials.cross_sections = '/Users/pablogarcia44/repo/cross_sections.xml'


#Geometry
def pincell(index_x,index_y):
    fuel_outer_radius = openmc.ZCylinder(r=0.39)
    clad_inner_radius = openmc.ZCylinder(r=0.40)
    clad_outer_radius = openmc.ZCylinder(r=0.46)
    pitch = 1.26
    left = openmc.XPlane(-pitch/2, boundary_type='transmission')
    right = openmc.XPlane(pitch/2, boundary_type='transmission')
    bottom = openmc.YPlane(-pitch/2, boundary_type='transmission')
    top = openmc.YPlane(pitch/2, boundary_type='transmission')
    fuel_region = -fuel_outer_radius
    gap_region = +fuel_outer_radius & -clad_inner_radius
    clad_region = +clad_inner_radius & -clad_outer_radius
    water_region = +left & -right & +bottom & -top & +clad_outer_radius
    fuel = openmc.Cell(name='fuel'+'_'+str(index_x)+'_'+str(index_y))
    fuel.fill = uo2
    fuel.region = fuel_region
    gap = openmc.Cell(name='gap'+'_'+str(index_x)+'_'+str(index_y))
    gap.region = gap_region
    clad = openmc.Cell(name='clad'+'_'+str(index_x)+'_'+str(index_y))
    clad.fill = zirconium
    clad.region = clad_region
    moderator = openmc.Cell(name='moderator'+'_'+str(index_x)+'_'+str(index_y))
    moderator.fill = water
    moderator.region = water_region 
    u = openmc.Universe(name='u'+'_'+str(index_x)+'_'+str(index_y),cells=(fuel, gap, clad, moderator))
    return(u,fuel,gap,clad,moderator)

def guide(index_x,index_y):
    clad_inner_radius = openmc.ZCylinder(r=0.5610)
    clad_outer_radius = openmc.ZCylinder(r=0.6020)
    pitch = 1.26
    left = openmc.XPlane(-pitch/2, boundary_type='transmission')
    right = openmc.XPlane(pitch/2, boundary_type='transmission')
    bottom = openmc.YPlane(-pitch/2, boundary_type='transmission')
    top = openmc.YPlane(pitch/2, boundary_type='transmission')       
    clad_region = +clad_inner_radius & -clad_outer_radius
    water_region = +left & -right & +bottom & -top & +clad_outer_radius
    water_guide_region= -clad_inner_radius
    water_guide = openmc.Cell(name='water_guide'+'_'+str(index_x)+'_'+str(index_y))
    water_guide.fill = water
    water_guide.region = water_guide_region
    clad = openmc.Cell(name='clad_guide'+'_'+str(index_x)+'_'+str(index_y))
    clad.fill = zirconium
    clad.region = clad_region
    moderator = openmc.Cell(name='moderator'+'_'+str(index_x)+'_'+str(index_y))
    moderator.fill = water
    moderator.region = water_region 
    u = openmc.Universe(name='u_guide'+'_'+str(index_x)+'_'+str(index_y),cells=(water_guide, clad, moderator))
    return(u,water_guide,water_guide,clad,moderator)
    
assembly = openmc.RectLattice()
pitch=1.26
dr=2e-1 # cm of water that is outside assembly
size=17 # #size of the assembly 
pitch_assembly=size*pitch+2*dr 
assembly.pitch = (pitch,pitch)
assembly.lower_left = (-size/2*pitch, -size/2*pitch)
#template 17x17 for GT
template_x = np.array([5, 8, 11, 3, 13, 2, 5, 8, 11, 14, 2, 5, 8,
                           11, 14, 2, 5, 8, 11, 14, 3, 13, 5, 8, 11])
template_y = np.array([2, 2, 2, 3, 3, 5, 5, 5, 5, 5, 8, 8, 8, 8,
                           8, 11, 11, 11, 11, 11, 13, 13, 14, 14, 14])
A= np.empty((size, size), dtype=openmc.universe.Universe)
for ix in range(size):
    for iy in range(size):
        indices_x = np.where(template_x == ix)[0]
        indices_y = np.where(template_y == iy)[0]
        if any(i in indices_y for i in indices_x):
            A[ix][iy]=guide(ix,iy)[0]
        else:
            A[ix][iy]=pincell(ix,iy)[0]      
assembly.universes = A
moderator_outside = openmc.Cell()
moderator_outside.fill = water
all_water = openmc.Universe()
all_water.add_cell(moderator_outside)
assembly.outer=all_water
min_x = openmc.XPlane(x0=-(size/2*pitch+dr), boundary_type='reflective')
max_x = openmc.XPlane(x0=0, boundary_type='reflective')
min_y = openmc.YPlane(y0=0, boundary_type='reflective')
max_y = openmc.YPlane(y0=+(size/2*pitch+dr), boundary_type='reflective')
min_z = openmc.ZPlane(z0=-10., boundary_type='reflective')
max_z = openmc.ZPlane(z0=+10., boundary_type='reflective')
root_cell = openmc.Cell(name='root cell', fill=assembly)
root_cell.region = +min_x & -max_x & +min_y & -max_y & +min_z & -max_z
model.geometry.root_universe = openmc.Universe(name='root universe')
model.geometry.root_universe.add_cell(root_cell)

# Settings
bounds = [-pitch_assembly/2, 0, -10, 0, +pitch_assembly/2, 10]
uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:], only_fissionable=True)
source = openmc.IndependentSource(space=uniform_dist)
source.time = openmc.stats.Uniform(0, 1e-6)
settings = openmc.Settings()
settings.source = source
settings = openmc.Settings()
settings.source = source
settings.batches = 500
settings.inactive = 100
settings.particles = 25000000
settings.output = {'tallies':True}
model.settings = settings

#Tallies
choose_group=openmc.mgxs.GROUP_STRUCTURES['XMAS-172']
name_group='XMAS172'
mgxs_lib = openmc.mgxs.Library(model.geometry)
groups = openmc.mgxs.EnergyGroups(choose_group)
mgxs_lib.energy_groups = groups
mgxs_lib.scatter_format = "legendre"
mgxs_lib.legendre_order = 3
mgxs_lib.mgxs_types = ['total', 'absorption', 'fission', 'nu-fission', 'chi', 'scatter matrix', 'nu-scatter matrix', 'inverse-velocity']
mgxs_lib.tally_trigger =openmc.Trigger('std_dev', 1e-2)
mgxs_lib.by_nuclide = False
mgxs_lib.domain_type = 'cell'
mgxs_lib.domains = model.geometry.get_all_material_cells().values()
mgxs_lib.build_library()
tallies = openmc.Tallies()
mgxs_lib.add_to_tallies_file(tallies, merge=True)
model.tallies = tallies

#model.export_to_xml()
model.export_to_model_xml()

#sp_file = model.run()