import mm_kinetic_analysis
import numpy as np
import math
import ast
from mm_kinetic_analysis import get_parser, Import_Kinetic_Data, MM_Kinetic_Solver, get_inputs, graph_kinetic_data

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

with open('substrate_data.txt', 'r') as file:
    lines = [line.strip() for line in file.readlines()]
    line_1 = int(lines[0])
    line_2 = float(lines[1])
    line_3 = float(lines[2])
    substrate_vals = [line_1, line_2, line_3]

inputs = get_inputs()
substrate = inputs.gen_substrate(substrate_vals)

with open('path_data.txt', 'r') as file:
    path = [line.strip() for line in file.readlines()]
    print(path)

data = Import_Kinetic_Data(path[0], substrate)

with open('column_data.txt', 'r') as file:
    content = file.read().strip()
    lines = ast.literal_eval(content)
    line_1 = int(lines[0])+2
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

sum_value_guess = []
sum_value_min = []
kinetic_parameters_all = []

with open('data_type.txt', 'r') as file:
    lines = [line.strip() for line in file.readlines()]
    line_1 = lines[0]
    data_type = [line_1]
    if len(lines[1]) > 0:
        line_2 = int(lines[1])
        data_type.append(line_2)


for i in range(len(vvalues_all)):
    vvalues = vvalues_all[i]
    if data_type[0] == "Absorbance":
        vvalues_abs = []
        for i in range(len(vvalues)):
            vv = vvalues[i]/(int(data_type[1]))
            vvalues_abs.append(vv)
        vvalues = vvalues_abs
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
        sum_value_min.append(s_min)
        kinetic_parameters_all.append(sol)
        print(f"Done Calculating Kinetic Parameters at V Value {i}")

best_v = np.min(sum_value_min)
sum_value_min = np.array(sum_value_min)
ind_min = np.where(sum_value_min == best_v)
ind_min = ind_min[0].astype(int)
kinetic_parameters = kinetic_parameters_all[ind_min[0]]
vvalues = vvalues_all[ind_min[0]]
if data_type[0] == "Absorbance":
    vvalues_abs = []
    for i in range(len(vvalues)):
        vv = vvalues[i]/(int(data_type[1]))
        vvalues_abs.append(vv)
    vvalues = vvalues_abs   


vval_calc = []
for i in range(len(substrate)):
    val = MM_Kinetic_Solver(kinetic_parameters[0], kinetic_parameters[1])
    calc = val.mm_equation(i, substrate)
    vval_calc.append(calc)

with open('name_data.txt', 'r') as file:
    name = [line.strip() for line in file.readlines()]
    print(name)

with open(f'{name[0]}.txt', 'w') as file:
    file.write(str(1))
    file.write("\n")
    file.write(str(kinetic_parameters[0]))
    file.write("\n")
    file.write(str(kinetic_parameters[1]))
    file.write("\n")



plot = graph_kinetic_data(name[0], substrate, vvalues, vval_calc, kinetic_parameters, 0)
plot.mm_graph()
