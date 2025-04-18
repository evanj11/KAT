#!/usr/bin/env python
# coding: utf-8

# In[1]:


import hill_kinetic_analysis
import numpy as np
import math
from hill_kinetic_analysis import get_parser, Import_Kinetic_Data, Hill_Kinetic_Solver, get_inputs, graph_kinetic_data


# In[2]:


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

with open('substrate_data.txt', 'r') as file:
    lines = [line.strip() for line in file.readlines()]
    print(lines)
    line_1 = int(lines[0])
    line_2 = float(lines[1])
    line_3 = int(lines[2])
    substrate_vals = [line_1, line_2, line_3]

inputs = get_inputs()
substrate = inputs.gen_substrate(substrate_vals)

with open('path_data.txt', 'r') as file:
    path = file.readlines()
    print(path)

data = Import_Kinetic_Data(path[0], substrate)

with open('column_data.txt', 'r') as file:
    lines = [line.strip() for line in file.readlines()]
    line_1 = int(lines[0])+3
    line_2 = int(lines[1])+3
    columns = [line_1, line_2]
    print(columns)

df = data.import_data(columns)

with open('time_data.txt', 'r') as file:
    lines = [line.strip() for line in file.readlines()]
    line_1 = int(lines[0])
    line_2 = int(lines[1])
    line_3 = int(lines[2])
    time = [line_1, line_2, line_3]

vvalues_all = data.gen_vvalues(df, time_min=time[0], time_max=time[1], steps=time[2])
print(vvalues_all)


# In[7]:


sum_value_guess = []
sum_value_min = []
vvalues = []
vv_std = []
for j in range(len(vvalues_all[0])):
    vvalue_it = []
    for i in range(len(vvalues_all)):
        vvalue_it.append(vvalues_all[i][j])
    vvalues.append(np.average(vvalue_it))
    vv_std.append(np.std(vvalue_it))

print(vv_std)
vm = (vvalues[0] + vvalues[1] + vvalues[2])/3
hv = vm/2
hv = int(hv)
vkm = find_nearest(vvalues, hv)
ind = np.where(vvalues == vkm)
ind = ind[0]
ind = ind.astype(int)
if ind.size == 0:
    ind = 0
    print("One or More V Values are NaN, Move on to Next V Value")
else:
    val = Hill_Kinetic_Solver(2, vm, substrate[ind[0]+1])
    s = val.sums(2, vm, substrate[ind[0]+1], vvalues, substrate)
    sum_value_guess.append(s)
    eq_to_min = val.full_equation(substrate, vvalues)
    df_dvmax, df_dh, df_dkm = val.partial_diff(eq_to_min)
    sol = val.minimize(df_dvmax, df_dh, df_dkm)
    val_min = Hill_Kinetic_Solver(sol[0], sol[1], sol[2])
    s_min = val_min.sums(sol[0], sol[1], sol[2], vvalues, substrate)
    kinetic_parameters = sol
    print(f"Done Calculating Kinetic Parameters")


# In[8]:


print(kinetic_parameters, s_min, vvalues)


# In[9]:


vval_calc = []
for i in range(len(substrate)):
    val = Hill_Kinetic_Solver(kinetic_parameters[0], kinetic_parameters[1], kinetic_parameters[2])
    calc = val.hill_equation(i, substrate)
    vval_calc.append(calc)


# In[10]:


print(vvalues[0])
spx, spy = inputs.linear_hill_xy(vvalues, substrate)


# In[13]:


poly1d_fn, linregx = inputs.linreg(spx, spy, 4, 30)


with open('entry_data.txt', 'r') as file:
    name = file.readlines()
    print(name)

plot = graph_kinetic_data(name[0], substrate, vvalues, vval_calc, kinetic_parameters, vv_std)
plot.with_inset(spx, spy, linregx, poly1d_fn)



