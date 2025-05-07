import mm_kinetic_analysis
import numpy as np
import math
from mm_kinetic_analysis import get_parser, Import_Kinetic_Data, MM_Kinetic_Solver, get_inputs, graph_kinetic_data


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
    line_1 = int(lines[0])+2
    line_2 = int(lines[1])+2
    columns = [line_1, line_2]
    print(columns)

df = data.import_data(columns)

vvalues_all = data.gen_vvalues(df, time_min=5, time_max=15, steps=10)
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
    val = MM_Kinetic_Solver(vm, substrate[ind[0]+1])
    s = val.sums(vm, substrate[ind[0]+1], vvalues, substrate)
    sum_value_guess.append(s)
    eq_to_min = val.full_equation(substrate, vvalues)
    df_dvmax, df_dkm = val.partial_diff(eq_to_min)
    sol = val.minimize(df_dvmax, df_dkm)
    val_min = MM_Kinetic_Solver(sol[0], sol[1])
    s_min = val_min.sums(sol[0], sol[1], vvalues, substrate)
    kinetic_parameters = sol
    print(f"Done Calculating Kinetic Parameters")


# In[8]:


print(kinetic_parameters, s_min, vvalues)


# In[9]:


vval_calc = []
for i in range(len(substrate)):
    val = MM_Kinetic_Solver(kinetic_parameters[0], kinetic_parameters[1])
    calc = val.mm_equation(i, substrate)
    vval_calc.append(calc)


with open('entry_data.txt', 'r') as file:
    name = file.readlines()
    print(name)

plot = graph_kinetic_data(name[0], substrate, vvalues, vval_calc, kinetic_parameters, vv_std)
plot.mm_graph()



