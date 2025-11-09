import sys
import mm_kinetic_analysis
import numpy as np
import math
import os
import pandas as pd
import ast
import argparse
import matplotlib.pyplot as plt
from scipy.optimize import minimize, differential_evolution
from mm_kinetic_analysis import get_parser, Import_Kinetic_Data, MM_Kinetic_Solver
from hill_kinetic_analysis import Hill_Kinetic_Solver
from mwc_kinetic_analysis import MWC_Kinetic_Solver, get_inputs, graph_kinetic_data
from knf_kinetic_analysis import KNF_Kinetic_Solver

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

def no_inset(model, name, substrate, vvalues, vval_calc, kinetic_parameters, vv_std):
    if model == "mm":
        fig, ax = plt.subplots(figsize=(4,3), dpi=250)
        plt.subplots_adjust(left=0.15, wspace=0.3, bottom=0.15)
        ax.errorbar(substrate, vvalues, yerr=vv_std, fmt="*", color='blue', label="Data", markersize=6, elinewidth=1, capsize=1.5, barsabove=True)
        ax.errorbar(substrate, vval_calc, fmt="o-", color='black', label="Calculated", markersize=2)
        ax.set_ylabel("V\u2080")
        ax.set_xlabel("[Substrate]")
        ax.set_title("M-M Kinetic Plot")
        plt.savefig(f'{name}.png')
        plt.savefig(f'{name}.svg')
    if model == "hill":
        fig, ax = plt.subplots(figsize=(4,3), dpi=250)
        plt.subplots_adjust(left=0.15, wspace=0.3, bottom=0.15)
        ax.errorbar(substrate, vvalues, yerr=vv_std, fmt="*", color='blue', label="Data", markersize=6, elinewidth=1, capsize=1.5, barsabove=True)
        ax.errorbar(substrate, vval_calc, fmt="o-", color='black', label="Calculated", markersize=2)
        ax.set_ylabel("V\u2080")
        ax.set_xlabel("[Substrate]")
        ax.set_title("Hill Kinetic Plot")
        ax.set_xscale("log")
        hill = '%.2f'%(kinetic_parameters[0])
        ax.annotate(f"Hill Coefficient = {hill}", xy=(5, 0), fontsize=7, fontstyle='italic')
        fig.savefig(f'{name}.png', dpi=300)
        fig.savefig(f'{name}.svg')
    if model == "mwc":
        fig, ax = plt.subplots(figsize=(4,3), dpi=250)
        plt.subplots_adjust(left=0.15, wspace=0.3, bottom=0.15)
        ax.errorbar(substrate, vvalues, yerr=vv_std, fmt="*", color='blue', label="Data", markersize=6, elinewidth=1, capsize=1.5, barsabove=True)
        ax.errorbar(substrate, vval_calc, fmt="o-", color='black', label="Calculated", markersize=2)
        ax.set_ylabel("V\u2080")
        ax.set_xlabel("[Substrate]")
        ax.set_title("MWC Kinetic Plot")
        ax.set_xscale("log")
        hill = '%.2f'%(kinetic_parameters[5])
        ax.annotate(f'# of Allosteric Sites = {hill}', 
            xy=(0.05, 0.97), # Adjust based on your bbox_to_anchor for the legend
            xycoords='axes fraction', 
            xytext=(0, 0), # Offset from xy in points
            textcoords='offset points',
            horizontalalignment='left', 
            verticalalignment='top',
            fontsize=7,
            fontstyle='italic')
        fig.savefig(f'{name}.png', dpi=300)
        fig.savefig(f'{name}.svg')
    if model == "knf":
        fig, ax = plt.subplots(figsize=(4,3), dpi=250)
        plt.subplots_adjust(left=0.15, wspace=0.3, bottom=0.15)
        ax.errorbar(substrate, vvalues, yerr=vv_std, fmt="*", color='blue', label="Data", markersize=6, elinewidth=1, capsize=1.5, barsabove=True)
        ax.errorbar(substrate, vval_calc, fmt="o-", color='black', label="Calculated", markersize=2)
        ax.set_ylabel("V\u2080")
        ax.set_xlabel("[Substrate]")
        ax.set_title("KNF Kinetic Plot")
        ax.set_xscale("log")
        hill = '%.2f'%(kinetic_parameters[3])
        ax.annotate(f'Gamma = {hill}', 
            xy=(0.05, 0.97), # Adjust based on your bbox_to_anchor for the legend
            xycoords='axes fraction', 
            xytext=(0, 0), # Offset from xy in points
            textcoords='offset points',
            horizontalalignment='left', 
            verticalalignment='top',
            fontsize=7,
            fontstyle='italic')
        fig.savefig(f'{name}.png', dpi=300)
        fig.savefig(f'{name}.svg')
    return fig

def parse_list_or_file(value):
    if os.path.isfile(value):
        with open(value) as f:
            lines = [line.strip() for line in f if line.strip()]
            return [float(x) for x in lines]
    else:
        return [float(x) for x in value.strip("[]").split(",")]
        
def main():
    parser = argparse.ArgumentParser(description="Welcome to KAT, please provide substrate conc., velocity values, model, and output filename")
    parser.add_argument("-s", "--substrate", type=str, nargs="+", required=True, help="List or file of substrate concentrations")
    parser.add_argument("-v", "--velocity", type=str, nargs="+", required=True, help="List or file containing pre-calculated velocity values")
    parser.add_argument("-m", "--model", required=True, choices=["mm", "hill", "mwc", "knf"], help="Choose which kinetic model to use")
    parser.add_argument("-o", "--output", type=str, default="kinetic_model_fitting", help="Output filename for graph and kinetic parameters")
    parser.add_argument("--calc-sub-dil", action="store_true", help=("Calculate substrate dilutions from substrate input file with format:\n"
                                                                     "1. Number of Substrate Dilutions\n"
                                                                     "2. Substrate Dilution Factor\n"
                                                                     "3. Highest Substrate Concentration\n"))
    parser.add_argument("--sites", type=int, help="Number of binding sites (required for KNF model only)")

    args = parser.parse_args()

    if args.model == "knf" and args.sites is None:
        parser.error("\033[91mThe --sites argument is required when using the KNF model\033[0m")
    
    print(f"\033[94mRunning {args.model} Kinetic Analysis\033[0m")

    if args.calc_sub_dil is True:
        sub_arg = args.substrate[0] if isinstance(args.substrate, list) else args.substrate
        if os.path.isfile(sub_arg):
            with open(sub_arg, 'r') as file:
                lines = [line.strip() for line in file.readlines()]
                line_1 = int(lines[0])
                line_2 = float(lines[1])
                line_3 = float(lines[2])
                substrate_vals = [line_1, line_2, line_3]
        else:
            substrate_vals = [float(x) for x in args.substrate]
    
        inputs = get_inputs()
        substrate = inputs.gen_substrate(substrate_vals)
    else:
        sub_arg = args.substrate[0] if isinstance(args.substrate, list) else args.substrate
        if os.path.isfile(sub_arg):
            with open(sub_arg, 'r') as file:
                lines = [line.strip() for line in file.readlines()]
                substrate = [float(line) for line in lines]
        else:
            substrate = [float(x) for x in args.substrate]

    print(f"Substrate Values are {substrate}")
    vvalues = []
    for v in args.velocity:
        vvalues.extend(parse_list_or_file(v))
    
    sum_value_guess = []
    sum_value_min = []
    if args.model == "mm":
        if vvalues[0] > vvalues[1]:
            vm = (vvalues[0] + vvalues[1] + vvalues[2]) / 3
        else:
            vm = (vvalues[-3:][0] + vvalues[-3:][1] + vvalues[-3:][2]) / 3
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
            val = MM_Kinetic_Solver(vm, substrate[ind[0]])
            s = val.sums(vm, substrate[ind[0]], vvalues, substrate)
            sum_value_guess.append(s)
            eq_to_min = val.full_equation(substrate, vvalues)
            df_dvmax, df_dkm = val.partial_diff(eq_to_min)
            sol = val.minimize(df_dvmax, df_dkm)
            val_min = MM_Kinetic_Solver(sol[0], sol[1])
            s_min = val_min.sums(sol[0], sol[1], vvalues, substrate)
            kinetic_parameters = sol
            print(f"Done Calculating Kinetic Parameters")
        vval_calc = []
        for i in range(len(substrate)):
            val = MM_Kinetic_Solver(kinetic_parameters[0], kinetic_parameters[1])
            calc = val.mm_equation(i, substrate)
            vval_calc.append(calc)
    elif args.model == "hill":
        if vvalues[0] > vvalues[1]:
            vm = (vvalues[0] + vvalues[1] + vvalues[2]) / 3
        else:
            vm = (vvalues[-3:][0] + vvalues[-3:][1] + vvalues[-3:][2]) / 3
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
            val = Hill_Kinetic_Solver(2, vm, substrate[ind[0]])
            s = val.sums(2, vm, substrate[ind[0]], vvalues, substrate)
            sum_value_guess.append(s)
            eq_to_min = val.full_equation(substrate, vvalues)
            df_dvmax, df_dh, df_dkm = val.partial_diff(eq_to_min)
            sol = val.minimize(df_dvmax, df_dh, df_dkm)
            val_min = Hill_Kinetic_Solver(sol[0], sol[1], sol[2])
            s_min = val_min.sums(sol[0], sol[1], sol[2], vvalues, substrate)
            kinetic_parameters = sol
            print(f"Done Calculating Kinetic Parameters", file=sys.stdout, flush=True)
        vval_calc = []
        for i in range(len(substrate)):
            val = Hill_Kinetic_Solver(kinetic_parameters[0], kinetic_parameters[1], kinetic_parameters[2])
            calc = val.hill_equation(i, substrate)
            vval_calc.append(calc)
    elif args.model == "mwc":
        V_R_guess = vvalues[0]
        V_T_guess = 0.05 * V_R_guess
        K_R_guess = substrate[np.abs(vvalues[0] - 0.5 * V_R_guess).argmin()]
        K_T_guess = 5 * K_R_guess   
        val = MWC_Kinetic_Solver(V_T_guess, V_R_guess, K_T_guess, K_R_guess, 100.0, 2)      
        p0 = [V_T_guess, V_R_guess, K_T_guess, K_R_guess, 100.0, 2.0]
        bounds = [(0, 1000), (0, 1000), (1e-6, 1000), (1e-6, 1000), (1e-3, 500), (0.5, 14)]
        result = differential_evolution(
                val.loss, 
                args=(substrate, vvalues),
                bounds=bounds,
                strategy='best1bin',
                maxiter=1000,
                popsize=20,
                tol=1e-7
            )
        refined = minimize(val.loss, result.x, args=(substrate, vvalues), bounds=bounds, method='L-BFGS-B', options={
            'maxiter': 10000,
            'ftol': 1e-12,    # tighter convergence on function value
            'gtol': 1e-10,    # tighter gradient convergence
            'eps': 1e-8       # step size
            })
        kinetic_parameters = refined.x
        s_min = refined.fun
        print(f"Done Calculating Kinetic Parameters", file=sys.stdout, flush=True)
        val_bf = MWC_Kinetic_Solver(kinetic_parameters[0], kinetic_parameters[1], kinetic_parameters[2], kinetic_parameters[3], kinetic_parameters[4], kinetic_parameters[5])
        vval_calc = val_bf.mwc_model(substrate)
    elif args.model == "knf":
        vmax_guess = vvalues[0]
        kd_guess = substrate[np.abs(vvalues[0] - 0.5 * vmax_guess).argmin()]
        kbasal_guess = np.min(vvalues)
        val = KNF_Kinetic_Solver(args.sites, vmax_guess, kd_guess, kbasal_guess, 2)
        p0 = [vmax_guess, kd_guess, kbasal_guess, 2.0]
        bounds = [(0, 10000), (0, 10000), (0, 10000), (-50, 50)] 
        result = differential_evolution(
                val.loss, 
                args=(substrate, vvalues),
                bounds=bounds,
                strategy='best1bin',
                maxiter=1000,
                popsize=20,
                tol=1e-7
                )
        refined = minimize(val.loss, result.x, args=(substrate, vvalues), bounds=bounds, method='L-BFGS-B', options={
            'maxiter': 10000,
            'ftol': 1e-12,    # tighter convergence on function value
            'gtol': 1e-10,    # tighter gradient convergence
            'eps': 1e-8       # step size
        }) 
        kinetic_parameters = refined.x
        s_min = refined.fun
        print(f"Done Calculating Kinetic Parameters", file=sys.stdout, flush=True)
        val_bf = KNF_Kinetic_Solver(args.sites, kinetic_parameters[0], kinetic_parameters[1], kinetic_parameters[2], kinetic_parameters[3])
        vval_calc = val_bf.knf_model(substrate)
    else:
        print("Please choose a model from [mm, hill, mwc, knf]")

    print(f'Kinetic Parameters are: {kinetic_parameters}\n', file=sys.stdout, flush=True)
    print(f'Minimum Residual Square Sum is {s_min}\n', file=sys.stdout, flush=True)
    vvalues_txt = [float(x) for x in vvalues]
    form_vv = ['%.4f' % elem for elem in vvalues_txt]
    print(f'Velocity values are {form_vv}\n', file=sys.stdout, flush=True)

    filename = args.output
    if not filename.endswith(".txt"):
        filename_txt = filename + ".txt"
    
    with open(filename_txt, 'w') as file:
        file.write(str(1))
        file.write("\n")
        for i in range(len(kinetic_parameters)):
            file.write(str(kinetic_parameters[i]))
            file.write("\n")

    plot = no_inset(args.model, args.output, substrate, vvalues, vval_calc, kinetic_parameters, 0)
    print("\033[94mScript finished successfully\033[0m")
    print("\033[1;35mThank you for using KAT!\033[0m")
    print("\033[1;35mPlease cite @evanj11 on GitHub\033[0m")

if __name__ == "__main__":
    main()


