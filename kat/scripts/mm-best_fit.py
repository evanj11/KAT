import kinetic_analysis
import numpy as np
import matplotlib.pyplot as plt
import math
from mm_kinetic_analysis import get_parser, Import_Kinetic_Data, MM_Kinetic_Solver, get_inputs, graph_kinetic_data

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

with open('substrate_data.txt', 'r') as file:
    lines = [line.strip() for line in file.readlines()]
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

sum_value_guess = []
sum_value_min = []
kinetic_parameters_all = []

for i in range(len(vvalues_all)):
    vvalues = vvalues_all[i]
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
        df_dvmax, df_dh, df_dkm = val.partial_diff(eq_to_min)
        sol = val.minimize(df_dvmax, df_dh, df_dkm)
        val_min = MM_Kinetic_Solver(sol[0], sol[1])
        s_min = val_min.sums(sol[0], sol[1], vvalues, substrate)
        sum_value_min.append(s_min)
        kinetic_parameters_all.append(sol)
        print(f"Done Calculating Kinetic Parameters at V Value {i}")

best_v = np.min(sum_value_min)
sum_value_min = np.array(sum_value_min)
ind_min = np.where(sum_value_min == best_v)
ind_min = ind_min[0].astype(int)
kinetic_parameters = kinetic_parameters_all[ind_min[0]]
vvalues = vvalues_all[ind_min[0]]

vval_calc = []
for i in range(len(substrate)):
    val = MM_Kinetic_Solver(kinetic_parameters[0], kinetic_parameters[1])
    calc = val.mm_equation(i, substrate)
    vval_calc.append(calc)

spx, spy = inputs.linear_hill_xy(vvalues, substrate)

poly1d_fn, linregx = inputs.linreg(spx, spy, 6, 17)

with open('entry_data.txt', 'r') as file:
    name = file.readlines()
    print(name)

plot = graph_kinetic_data(name[0], substrate, vvalues, vval_calc, kinetic_parameters, 0)
plot.mm_graph()
