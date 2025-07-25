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
        line_3 = float(lines[2])
        substrate_vals = [line_1, line_2, line_3]

    inputs = get_inputs()
    substrate = inputs.gen_substrate(substrate_vals)

    with open(os.path.join(work_dir, 'path_rep_data.txt'), 'r') as file:
        pathds = [line.strip() for line in file.readlines()]
        print(pathds)
    with open(os.path.join(work_dir, 'blank_rep_data.txt'), 'r') as file:
        pathbs = [line.strip() for line in file.readlines()]
    
    col_max = len(substrate) + 2
    columns = [2, col_max]

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
    
    vvalues_rep = []
    kinetic_parameters_rep = []
    vval_calc_rep = []
    for i in range(len(pathds)):
        data = Import_Kinetic_Data(pathds[i], substrate)
        df_data = data.import_data(columns)
        if len(pathds) == len(pathbs):
            blank = Import_Kinetic_Data(pathbs[i], substrate)
            df_blank = blank.import_data(columns)
            df = pd.DataFrame(df_data.values - df_blank.values)
            print(df)
        elif 'None' not in pathbs and 1 == len(pathbs):
            blank = Import_Kinetic_Data(pathbs[0], substrate)
            df_blank = blank.import_data(columns)
            df = pd.DataFrame(df_data.values - df_blank.values)
        else:
            df = df_data
   
        if 'True' in time:
            lin_range = data.gen_lin_range(df, time[1])
            print(lin_range)
            step = lin_range[1] - lin_range[0]
            time = [lin_range[0], lin_range[1], step, line_2]

    
        vvalues_all = data.gen_vvalues(df, time_min=time[0], time_max=time[1], steps=time[2], v_win=time[3])
        print(vvalues_all)

        sum_value_guess = []
        sum_value_min = []
        vvalues = []
        vv_std = []
        for j in range(len(vvalues_all[0])):
            vvalue_it = []
            for av in range(len(vvalues_all)):
                vvalue_it.append(vvalues_all[av][j])
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
            print("One or More V Values are NaN, Move on to Next V Value")
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
            print(f"Done Calculating Kinetic Parameters")

        print(kinetic_parameters, s_min, vvalues)

        vval_calc = []
        for l in range(len(substrate)):
            val = Hill_Kinetic_Solver(kinetic_parameters[0], kinetic_parameters[1], kinetic_parameters[2])
            calc = val.hill_equation(l, substrate)
            vval_calc.append(calc)

        spx, spy = inputs.linear_hill_xy(vvalues, substrate)
        print(spx, spy)

        poly1d_fn, linregx = inputs.linreg(spx, spy)

        with open(os.path.join(work_dir, 'name_data.txt'), 'r') as file:
            name = [line.strip() for line in file.readlines()]
            print(name)

        with open(os.path.join(work_dir, f'{name[0]}_{i}.txt'), 'w') as file:
            file.write(str(kinetic_parameters[0]))
            file.write("\n")
            file.write(str(kinetic_parameters[1]))
            file.write("\n")
            file.write(str(kinetic_parameters[2]))
            file.write("\n")
        
        vvalues_rep.append(vvalues)
        vval_calc_rep.append(vval_calc)
        kinetic_parameters_rep.append(kinetic_parameters)
        plot = graph_kinetic_data(os.path.join(work_dir, f"{name[0]}_{i}"), substrate, vvalues, vval_calc, kinetic_parameters, 0)
        plot.with_inset(spx, spy, linregx, poly1d_fn)
    vval_sub = []
    vval_calc_sub = []
    for j in range(len(substrate)):
        vv_sub_temp = []
        vv_calc_temp = []
        vv_sub_temp.extend(float(vval[j]) for vval in vvalues_rep)
        vv_calc_temp.extend(float(vv_calc[j]) for vv_calc in vval_calc_rep)
        vval_sub.append(vv_sub_temp)
        vval_calc_sub.append(vv_calc_temp)
    kinetic_parameters_sub = []
    print(kinetic_parameters_rep)
    for m in range(3):
        kp_temp = []
        kp_temp.extend(float(kp_sub[m]) for kp_sub in kinetic_parameters_rep)
        kinetic_parameters_sub.append(kp_temp)
    print(kinetic_parameters_sub)
    vval_rep_avg = []
    vval_rep_std = []
    vval_calc_rep_avg = []
    vval_calc_rep_std = []
    for vv_sub in vval_sub:
        avg = np.average(vv_sub)
        vval_rep_avg.append(avg)
        std = np.std(vv_sub)
        vval_rep_std.append(std)
    for vv_calc in vval_calc_sub:
        avg_calc = np.average(vv_calc)
        vval_calc_rep_avg.append(avg_calc)
        std_calc = np.std(vv_calc)
        vval_calc_rep_std.append(std_calc)
    
    kinetic_parameters_avg = []
    kinetic_parameters_std = []
    for kp in kinetic_parameters_sub:
        avg_kp = np.average(kp)
        kinetic_parameters_avg.append(avg_kp)
        std_kp = np.std(kp)
        kinetic_parameters_std.append(std_kp)
    print(kinetic_parameters_avg, kinetic_parameters_std)

    with open(os.path.join(work_dir, f'{name[0]}.txt'), 'w') as file:
        hill_avg = '%.3f'%(kinetic_parameters_avg[0])
        hill_std = '%.3f'%(kinetic_parameters_std[0])
        file.write(str(hill_avg))
        file.write(" \u00B1 ")
        file.write(str(hill_std))
        file.write("\n")
        vmax_avg = '%.3f'%(kinetic_parameters_avg[1])
        vmax_std = '%.3f'%(kinetic_parameters_std[1])
        file.write(str(vmax_avg))
        file.write(" \u00B1 ")
        file.write(str(vmax_std))
        file.write("\n")
        km_avg = '%.3f'%(kinetic_parameters_avg[2])
        km_std = '%.3f'%(kinetic_parameters_std[2])
        file.write(str(km_avg))
        file.write(" \u00B1 ")
        file.write(str(km_std))
        file.write("\n")

    hill_tot = hill_avg + '\u00B1' + hill_std
    plot_rep = graph_kinetic_data(os.path.join(work_dir, name[0]), substrate, vvalues_rep, vval_calc_rep, hill_tot, 0)
    with open(os.path.join(work_dir, "mutant.txt"), 'r') as file:
        mut = [line.strip() for line in file.readlines()]

    if 'mutant' in mut:
        with open(os.path.join(work_dir, f'{name[0]}_mutant_data.txt'), 'w') as file:
            for mut in range(len(kinetic_parameters_rep)):
                file.write(str("%.3f"%(kinetic_parameters_rep[mut][0])))
                file.write(",")
                file.write(str("%.3f"%(kinetic_parameters_rep[mut][1])))
                file.write(",")
                file.write(str("%.3f"%(kinetic_parameters_rep[mut][2])))
                file.write("\n")
        plot_rep.mut_rep(kinetic_parameters_rep)
    else:
        plot_rep.rep_no_inset(vval_rep_avg, vval_calc_rep_avg, vval_rep_std, vval_calc_rep_std)

if __name__ == "__main__":
    main()

