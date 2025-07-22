import hill_kinetic_analysis
import numpy as np
import pandas as pd
import math
import ast
import os
from hill_kinetic_analysis import get_parser, Import_Kinetic_Data, Hill_Kinetic_Solver, get_inputs, graph_kinetic_data


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

def main():
    work_dir = os.environ.get('WORKING_DIR')

    with open(os.path.join(work_dir, 'substrate_data.txt'), 'r') as file:
        lines = [line.strip() for line in file.readlines()]
        line_1 = int(lines[0])
        line_2 = float(lines[1])
        line_3 = int(lines[2])
        substrate_vals = [line_1, line_2, line_3]

    inputs = get_inputs()
    substrate = inputs.gen_substrate(substrate_vals)

    with open(os.path.join(work_dir, 'path_data.txt'), 'r') as file:
        path = [line.strip() for line in file.readlines()]
        print(path)
    with open(os.path.join(work_dir, 'column_data.txt'), 'r') as file:
        content = file.read().strip()  # Read entire content and remove whitespace
        lines = ast.literal_eval(content)
        line_1 = int(lines[0])+2
        line_2 = int(lines[1])+3
        columns = [line_1, line_2]
        print(columns)

    data = Import_Kinetic_Data(path[0], substrate)
    df_data = data.import_data(columns)
    if path[1].strip():
        blank = Import_Kinetic_Data(path[1], substrate)
        df_blank = blank.import_data(columns)
        df = pd.DataFrame(df_data.values - df_blank.values)
    else:
        df = df_data
    with open(os.path.join(work_dir, 'time_data.txt'), 'r') as file:
        lines = [line.strip() for line in file.readlines()]
        line_1 = str(lines[0])
        line_2 = int(lines[1])
        if len(lines) > 2:
            line_3 = int(lines[2])
            line_4 = int(lines[3])
            line_5 = int(lines[4])
            time = [line_3, line_4, line_5, line_2]
        else:
            time = [line_1, line_2]
    
    if 'True' in time:
        lin_range = data.gen_lin_range(df, time[1])
        print(lin_range)
        step = lin_range[1] - lin_range[0]
        time = [lin_range[0], lin_range[1], step, line_2]

    vvalues_all = data.gen_vvalues(df, time_min=time[0], time_max=time[1], steps=time[2], v_win=time[3])
    print(vvalues_all)

    sum_value_guess = []
    sum_value_min = []
    kinetic_parameters_all = []
    with open(os.path.join(work_dir, 'data_type.txt'), 'r') as file:
        lines = [line.strip() for line in file.readlines()]
        line_1 = lines[0]
        data_type = [line_1]
        if 2 == len(lines):
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
            val = Hill_Kinetic_Solver(2, vm, substrate[ind[0]+1])
            s = val.sums(2, vm, substrate[ind[0]+1], vvalues, substrate)
            sum_value_guess.append(s)
            eq_to_min = val.full_equation(substrate, vvalues)
            df_dvmax, df_dh, df_dkm = val.partial_diff(eq_to_min)
            sol = val.minimize(df_dvmax, df_dh, df_dkm)
            val_min = Hill_Kinetic_Solver(sol[0], sol[1], sol[2])
            s_min = val_min.sums(sol[0], sol[1], sol[2], vvalues, substrate)
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
        val = Hill_Kinetic_Solver(kinetic_parameters[0], kinetic_parameters[1], kinetic_parameters[2])
        calc = val.hill_equation(i, substrate)
        vval_calc.append(calc)

    spx, spy = inputs.linear_hill_xy(vvalues, substrate)

    poly1d_fn, linregx = inputs.linreg(spx, spy)

    with open(os.path.join(work_dir, 'name_data.txt'), 'r') as file:
        name = [line.strip() for line in file.readlines()]
        print(name)

    with open(os.path.join(work_dir, f'{name[0]}.txt'), 'w') as file:
        file.write(str(kinetic_parameters[0]))
        file.write("\n")
        file.write(str(kinetic_parameters[1]))
        file.write("\n")
        file.write(str(kinetic_parameters[2]))
        file.write("\n")

    plot = graph_kinetic_data(os.path.join(work_dir, name[0]), substrate, vvalues, vval_calc, kinetic_parameters, 0)
    plot.with_inset(spx, spy, linregx, poly1d_fn)

if __name__ == '__main__':
    main()
