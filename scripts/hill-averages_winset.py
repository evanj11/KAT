import sys
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
    
    print("Running Hill Kinetic Analysis\n")
    with open(os.path.join(work_dir, 'substrate_data.txt'), 'r') as file:
        lines = [line.strip() for line in file.readlines()]
        line_1 = int(lines[0])
        line_2 = float(lines[1])
        line_3 = float(lines[2])
        substrate_vals = [line_1, line_2, line_3]

    inputs = get_inputs()
    substrate = inputs.gen_substrate(substrate_vals)

    with open(os.path.join(work_dir, 'path_data.txt'), 'r') as file:
        path = [line.strip() for line in file.readlines()]
    col_max = len(substrate) + 2
    columns = [2, col_max]

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
        print(f'Linear Range: {lin_range}\n', file=sys.stdout, flush=True)
        step = lin_range[1] - lin_range[0]
        time = [lin_range[0], lin_range[1], step, line_2]

    
    vvalues_all = data.gen_vvalues(df, time_min=time[0], time_max=time[1], steps=time[2], v_win=time[3])

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

    with open(os.path.join(work_dir, 'data_type.txt'), 'r') as file:
        lines = [line.strip() for line in file.readlines()]
        line_1 = lines[0]
        data_type = [line_1]
        if 2 == len(lines):
            line_2 = int(lines[1])
            data_type.append(line_2)

    if data_type[0] == "Absorbance":
        vvalues_abs = []
        for i in range(len(vvalues)):
            vv = vvalues[i] / (int(data_type[1]))
            vvalues_abs.append(vv)
        vvalues = vvalues_abs

    vm = (vvalues[0] + vvalues[1] + vvalues[2]) / 3
    hv = vm / 2
    hv = int(hv)
    vkm = find_nearest(vvalues, hv)
    ind = np.where(vvalues == vkm)
    ind = ind[0]
    ind = ind.astype(int)
    if ind.size == 0:
        ind = 0
        print("One or More V Values are NaN, Move on to Next V Value", file=sys.stdout, flush=True)
    else:
        val = Hill_Kinetic_Solver(2, vm, substrate[ind[0] + 1])
        s = val.sums(2, vm, substrate[ind[0] + 1], vvalues, substrate)
        sum_value_guess.append(s)
        eq_to_min = val.full_equation(substrate, vvalues)
        df_dvmax, df_dh, df_dkm = val.partial_diff(eq_to_min)
        sol = val.minimize(df_dvmax, df_dh, df_dkm)
        val_min = Hill_Kinetic_Solver(sol[0], sol[1], sol[2])
        s_min = val_min.sums(sol[0], sol[1], sol[2], vvalues, substrate)
        kinetic_parameters = sol
        print(f"Done Calculating Kinetic Parameters", file=sys.stdout, flush=True)

    print(f'Kinetic Parameters are {kinetic_parameters}\n', file=sys.stdout, flush=True)
    print(f'Minimum Residual Sum is {s_min}\n', file=sys.stdout, flush=True)
    vvalues_txt = [float(x) for x in vvalues]
    form_vv = ['%.4f' % elem for elem in vvalues_txt]
    print(f'Velocity values are {form_vv}\n', file=sys.stdout, flush=True)

    vval_calc = []
    for i in range(len(substrate)):
        val = Hill_Kinetic_Solver(kinetic_parameters[0], kinetic_parameters[1], kinetic_parameters[2])
        calc = val.hill_equation(i, substrate)
        vval_calc.append(calc)

    spx, spy = inputs.linear_hill_xy(vvalues, substrate)

    poly1d_fn, linregx = inputs.linreg(spx, spy)

    with open(os.path.join(work_dir, 'name_data.txt'), 'r') as file:
        name = [line.strip() for line in file.readlines()]

    with open(os.path.join(work_dir, f'{name[0]}.txt'), 'w') as file:
        file.write(str(kinetic_parameters[0]))
        file.write("\n")
        file.write(str(kinetic_parameters[1]))
        file.write("\n")
        file.write(str(kinetic_parameters[2]))
        file.write("\n")

    plot = graph_kinetic_data(os.path.join(work_dir, name[0]), substrate, vvalues, vval_calc, kinetic_parameters, vv_std)
    plot.with_inset(spx, spy, linregx, poly1d_fn)


if __name__ == "__main__":
    main()

