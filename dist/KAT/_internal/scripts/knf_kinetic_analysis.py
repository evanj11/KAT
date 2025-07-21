import math
import scipy as sc
import numpy as np
import pandas as pd
import argparse
from math import comb
from scipy.optimize import fsolve, least_squares, Bounds, minimize
from sympy import symbols, diff, solve, nsolve, checksol
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

def get_parser():
    """"Return a command line parser for this script."""
    parser = argparse.ArgumentParser(
        description="This script reads a .csv file of kinetic data (fluorescence) "
        "and returns a graph of the data along with Hill coefficient")
    parser.add_argument(
        "-f",
        "--file-name",
        dest="file_name",
        required=True,
        help="CSV file containing fluorescent data "
        )
    parser.add_argument(
        "-w",
        "--wells",
        dest="wells",
        required=True,
        type=int,
        nargs="+",
        help="Where is the data located? Please provide comma deliminated list "
        "(4-13) for Wells A, (16-25) for Wells B, (28-37) for Wells C "
        )
    parser.add_argument(
        "-s",
        "--substrate",
        dest="substrate",
        type=float,
        nargs="+",
        required=True,
        help="List of #Substrate Concetrations tested, Dilution Factor used, Max Concentration "
        )
    parser.add_argument(
        "-o",
        "--output",
        dest="output",
        required=True,
        default="hillplot",
        help="Output filename for graphs "
        )

    return parser

def find_first_and_last_small_diff(data, threshold=10.0):
    first_index = None
    last_index = None

    for i in range(len(data) - 1): 
        x1 = data[i]
        x2 = data[i + 1]

        if x1 == 0:
            continue  # Avoid division by zero

        percent_diff = abs((x2 - x1) / x1) * 100 

        if percent_diff < threshold:
            if first_index is None:
                first_index = i + 1  # First match
            last_index = i + 1      # Update last match

    if first_index is not None and last_index is not None:
        return [first_index, last_index]
    else:
        return None, None  # No matches found

class KNF_Kinetic_Solver:
    """
    Kinetic_Solver( hill coefficient, v max value, km value)
    """

    def __init__(self, sites, vmax, kd, k_basal, gamma):
        self.sites = sites
        self.vmax = vmax
        self.kd = kd
        self.k_basal = k_basal
        self.gamma = gamma
        
    def knf_model(self, S):
        vmax = self.vmax
        kd = self.kd
        k_basal = self.k_basal
        gamma = self.gamma
        self.S = S
        sites = int(self.sites)

        E_total = 1.0
        numerator = 0.0
        denominator = 0.0
        
        for i in range(sites+1):
            C_i = comb(sites, i) * (S / kd)**i
            if i == 0:
                k_i = 0
            else:
                k_i = k_basal + (vmax - k_basal) * (i / sites)**gamma
            numerator += k_i * C_i
            denominator += C_i
        return E_total * numerator / denominator

    def loss(self, params, S_data, v_data):
        sites = int(self.sites)
        vmax, kd, k_basal, gamma = params
        ks = KNF_Kinetic_Solver(sites, vmax, kd, k_basal, gamma)
        v_model = ks.knf_model(S_data)
        residual = np.sum((v_data - v_model) ** 2)
        if gamma > 30:
            penalty_y = 10 * (gamma - 8) ** 2
        else:
            penalty_y = 0

        return residual + penalty_y

    def square_sum(self, cal, dat):
        x = cal
        y = dat
        eq2 = (x-y)**2
        return eq2

    def sums(self, V_T, V_R, K_T, K_R, L0, n, vvalues, substrate):
        self.V_T = V_T
        self.V_R = V_R
        self.K_T = K_T
        self.K_R = K_R
        self.L0 = L0
        self.n = n
        self.substrate = substrate
        self.vvalues = vvalues
        sums = []
        for i in range(len(vvalues)):
            ks = MWC_Kinetic_Solver(V_T, V_R, K_T, K_R, L0, n)
            s = ks.square_sum(vvalues[i], ks.mwc_equation(i, substrate))
            sums.append(s)
        value = sum(sums)
        return value

    def full_equation(self, substrate, vvalues):
        equation = 0
        self.substrate = substrate
        self.vvalues = vvalues 
        for i in range(len(substrate)):
            V_T, V_R, K_T, K_R, L0, n = symbols('V_T V_R K_T K_R L0 n')
            sub = substrate[i]
            vval = vvalues[i]
            term_T = (1 + sub / K_T) ** n
            term_R = (1 + sub / K_R) ** n
            denom = L0 * term_T + term_R
            eq1 = (vval - ((V_T * L0 * term_T + V_R * term_R) / denom)) ** 2
            equation += eq1
        return equation

    def partial_diff(self, equation):
        V_T, V_R, K_T, K_R, L0, n = symbols('V_T V_R K_T K_R L0 n')
        f = equation
        df_dV_T = diff(f, V_T)
        df_dV_R = diff(f, V_R)
        df_dK_T = diff(f, K_T)
        df_dK_R = diff(f, K_R)
        df_dL0 = diff(f, L0)
        df_dn = diff(f, n)
        return df_dV_T, df_dV_R, df_dK_T, df_dK_R, df_dL0, df_dn

    def minimize(self, df_dV_T, df_dV_R, df_dK_T, df_dK_R, df_dL0, df_dn):
        V_T, V_R, K_T, K_R, L0, n = symbols('V_T V_R K_T K_R L0 n')
        eq1 = df_dV_T
        eq2 = df_dV_R
        eq3 = df_dK_T
        eq4 = df_dK_R
        eq5 = df_dL0
        eq6 = df_dn
        sol = nsolve((eq1, eq2, eq3, eq4, eq5, eq6), (V_T, V_R, K_T, K_R, L0, n), (self.V_T, self.V_R, self.K_T, self.K_R, self.L0, self.n), prec=15, solver="bisect", verify=False)
        return [sol[0], sol[1], sol[2], sol[3], sol[4], sol[5]]

class Import_Kinetic_Data:
    """
    Import .csv file and compute v values
    """
    def __init__(self, fname, substrate):
        self.fname = fname
        self.substrate = substrate
    
    def import_data(self, columns):
        self.columns = columns
        col_min = int(columns[0])
        col_max = int(columns[1])
        fname = self.fname
        df = pd.read_csv(fname, encoding='ISO-8859-1', usecols=range(col_min, col_max), nrows=31)
        df = df.apply(pd.to_numeric, errors="coerce")
        return df

    def gen_vvalues(self, df, time_min=15, time_max=40, steps=25, v_win=5):
        self.time_min = time_min
        self.time_max = time_max
        self.steps = steps
        self.v_win = v_win
        substrate = self.substrate
        arr = df.to_numpy()
        vvalues = []
        for i in np.linspace(time_min, time_max, steps):
            vval_time = []
            for j in range(len(substrate)):
                i = int(i)
                k = int(i+v_win)
                v = abs(arr[i][j] - arr[k][j])/v_win
                vval_time.append(v)
            vvalues.append(vval_time)
        return vvalues
    
    def gen_lin_range(self, df, v_win):
        self.v_win = v_win
        arr = df.to_numpy()
        substrate = self.substrate
        lin_range_all = []
        velocities_all = []
        for i in range(len(arr) - v_win):
            velocities = []
            for j in range(len(substrate)):
                i = int(i)
                k = int(i+v_win)
                v = abs(arr[i][j] - arr[k][j])/v_win
                velocities.append(v)
            velocities_all.append(velocities)
        for k in range(len(substrate)):
            vals = [val[k] for val in velocities_all if None not in val]
            if k > 4:
                lin_range_sub = find_first_and_last_small_diff(vals, threshold=7.5)
            else:
                lin_range_sub = find_first_and_last_small_diff(vals, threshold=2.5)
            lin_range_all.append(lin_range_sub)
        print(lin_range_sub)
        print(lin_range_all)
        first_indices = [pair[0] for pair in lin_range_all if None not in pair]
        last_indices  = [pair[1] for pair in lin_range_all if None not in pair]
        avg_first = sum(first_indices) / len(first_indices)
        closest_first = min(first_indices, key=lambda x: abs(x - avg_first))
        avg_last = sum(last_indices) / len(last_indices)
        closest_last = min(last_indices, key=lambda x: abs(x - avg_last))
        return [closest_first, closest_last]

class get_inputs:
    def find_nearest(self, array, value):
        self.array = array
        self.value = value
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        return array[idx]

    def gen_substrate(self, sub):
        self.sub = sub
        sub_num = int(sub[0])
        dilution = sub[1]
        initial_sub = int(sub[2])
        substrate = [initial_sub]
        for i in range(1, sub_num):
            i = int(i)
            a = substrate[i-1]/dilution
            substrate.append(a)
        return substrate

    def linear_hill_xy(self, vvalues, substrate):
        self.vvalues = vvalues
        self.substrate = substrate
        spy = []
        spx = []
        for i in range(1, len(substrate)):
            if vvalues[i] == 0:
                print('skip')
            else:
                yval = (vvalues[i])/(vvalues[0]-vvalues[i])
            if yval <= 0:
                print('skip')
                print(f'yval at x=vvalues[{i}] < 0')
            else:
                x = math.log(substrate[i])
                spx.append(x)
                y = math.log(yval)
                spy.append(y)
        return spx, spy

    def linreg(self, spx, spy):
        self.spx = spx
        self.spy = spy
        lin_min = 2
        lin_max = len(spx) - 2 
        linregx = []
        linregy = []
        for i in range(lin_min,lin_max):
            linregx.append(spx[i])
            linregy.append(spy[i])
        poly1d_fn = np.poly1d(np.polyfit(linregx, linregy, 1))
        return poly1d_fn, linregx

class graph_kinetic_data:
    def __init__(self, name, substrate, vvalues, vval_calc, kinetic_parameters, vv_std):
        self.name = name
        self.substrate = substrate
        self.vvalues = vvalues
        self.vval_calc = vval_calc
        self.kinetic_parameters = kinetic_parameters
        self.vv_std = vv_std

    def with_inset(self, spx, spy, linregx, poly1d_fn, bbox=(220, 150, 600, 500), ax_fs=5):
        substrate = self.substrate
        vvalues = self.vvalues
        vval_calc = self.vval_calc
        kinetic_parameters = self.kinetic_parameters
        self.spx = spx
        self.spy = spy
        self.linregx = linregx
        self.poly1d_fn = poly1d_fn
        self.bbox = bbox
        self.ax_fs = ax_fs
        vv_std = self.vv_std

        name = self.name
        fig, ax = plt.subplots(figsize=(4,3), dpi=250)
        plt.subplots_adjust(left=0.15, wspace=0.3, bottom=0.15)
        ax.errorbar(substrate, vvalues, yerr=vv_std, fmt="*", color='blue', label="Data", markersize=6, elinewidth=1, capsize=1.5, barsabove=True)
        ax.errorbar(substrate, vval_calc, fmt="o-", color='black', label="Calculated", markersize=2)
        ax.set_ylabel("V\u2080")
        ax.set_xlabel("[Substrate]")
        ax.set_title("MWC Kinetic Plot")
        ax.set_xscale("log")
        fig.savefig(f'{name}.png', dpi=300)
        fig.savefig(f'{name}.svg')
        return fig

    def no_inset(self):
        substrate = self.substrate
        vvalues = self.vvalues
        vval_calc = self.vval_calc
        kinetic_parameters = self.kinetic_parameters
        vv_std = self.vv_std
        
        name = self.name
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

    def rep_no_inset(self, vval_rep_avg, vval_calc_rep_avg, vval_rep_std, vval_calc_rep_std):
        substrate = self.substrate
        vvalues_rep = self.vvalues
        vval_calc_rep = self.vval_calc
        hill = self.kinetic_parameters
        vv_std = np.zeros_like(substrate)
        self.vval_rep_avg = vval_rep_avg
        self.vval_rep_std = vval_rep_std
        self.vval_calc_rep_avg = vval_calc_rep_avg
        self.vval_calc_rep_std = vval_calc_rep_std

        name = self.name
        colors = ['red', 'green', 'purple', 'pink', 'gray']
        labels = ['Rep1', 'Rep2', 'Rep3', 'Rep4', 'Rep5']
        fig, ax = plt.subplots(figsize=(4,3), dpi=250)
        plt.subplots_adjust(left=0.15, wspace=0.4, bottom=0.15)
        for i in range(len(vvalues_rep)):
            ax.errorbar(substrate, vvalues_rep[i], yerr=vv_std, fmt="*", color=colors[i], label=labels[i], markersize=6, elinewidth=1, capsize=1.5, barsabove=True, alpha=0.6)
        ax.errorbar(substrate, vval_rep_avg, yerr=vval_rep_std, fmt="d", color='black', label="Replicas Average", markersize=4, elinewidth=1, capsize=1.5, barsabove=True)
        ax.errorbar(substrate, vval_calc_rep_avg, yerr=vval_calc_rep_std, fmt="-", color='slateblue', label="Calculated Average", markersize=5, elinewidth=1, capsize=1.5, barsabove=True)
        ax.set_ylabel("V\u2080")
        ax.set_xlabel("[Substrate]")
        ax.set_title("KNF Kinetic Plot")
        ax.set_xscale("log")
        fig.legend(loc='upper left', bbox_to_anchor=(0.01, 0.94), bbox_transform=plt.gca().transAxes, fontsize=8)
        ax.annotate(f'Gamma = {hill}', 
            xy=(0.05, 0.97), # Adjust based on your bbox_to_anchor for the legend
            xycoords='axes fraction', 
            xytext=(0, 0), # Offset from xy in points
            textcoords='offset points',
            horizontalalignment='left', 
            verticalalignment='top',
            fontsize=7,
            fontstyle='italic')
        try:
            fig.savefig(f"{name}.png", dpi=300)
            fig.savefig(f"{name}.svg")
            print("Plot saved successfully.")
        except Exception as e:
            print("Error saving plot:", e)
        return fig
    
    def mut_rep(self, kinetic_parameters_mut):
        substrate = self.substrate
        vvalues_rep = self.vvalues
        vval_calc_rep = self.vval_calc
        self.kinetic_parameters_mut = kinetic_parameters_mut
        vv_std = np.zeros_like(substrate)

        name = self.name
        colors = ['red', 'green', 'mediumpurple', 'pink', 'lightsteelblue']
        colors_calc = ['darkred', 'darkgreen', 'orchid', 'pinkvioletred', 'slategrey']
        labels = ['Rep1', 'Rep2', 'Rep3', 'Rep4', 'Rep5']
        fig, ax = plt.subplots(figsize=(4,3), dpi=250)
        plt.subplots_adjust(left=0.15, wspace=0.4, bottom=0.15)
        for i in range(len(vvalues_rep)):
            ax.errorbar(substrate, vvalues_rep[i], yerr=vv_std, fmt="*", color=colors[i], label=labels[i], markersize=6, elinewidth=1, capsize=1.5, barsabove=True, alpha=1.0)
            ax.errorbar(substrate, vval_calc_rep[i], fmt="o-", color=colors_calc[i], label="Calculated", markersize=2, alpha=1.0)
            hill = "%.3f"%(kinetic_parameters_mut[i][3])
            ax.annotate(f'Gamma = {hill}', 
            xy=(0.05, 0.97-0.05*i), # Adjust based on your bbox_to_anchor for the legend
            xycoords='axes fraction', 
            xytext=(0, 0), # Offset from xy in points
            textcoords='offset points',
            horizontalalignment='left', 
            verticalalignment='top',
            fontsize=7,
            fontstyle='italic')
  
        ax.set_ylabel("V\u2080")
        ax.set_xlabel("[Substrate]")
        ax.set_title("KNF Kinetic Plot")
        ax.set_xscale("log")
        scale = len(vvalues_rep)
        fig.legend(loc='upper left', bbox_to_anchor=(0.01, 0.94-0.04*scale), bbox_transform=plt.gca().transAxes, fontsize=8)
        try:
            fig.savefig(f"{name}.png", dpi=300)
            fig.savefig(f"{name}.svg")
            print("Plot saved successfully.")
        except Exception as e:
            print("Error saving plot:", e)
        return fig

    def lineweaver_burk(self):
        substrate = self.substrate
        vvalues = self.vvalues
        vval_calc = self.vval_calc
        vv_std = self.vv_std
        
        lb_sub = []
        lb_vv = []
        lb_vvc = []
        lb_vv_std = []
        for i in range(len(substrate)):
            sub = 1/substrate[i]
            lb_sub.append(sub)
            vv = 1/vvalues[i]
            lb_vv.append(vv)
            cal = 1/vval_calc[i]
            lb_vvc.append(cal)
            std = 1/vv_std[0]
            lb_vv_std.append(std)

        name = self.name
        fig, ax = plt.subplots(figsize=(4,3), dpi=250)
        plt.subplots_adjust(left=0.15, wspace=0.3, bottom=0.15)
        ax.errorbar(lb_sub, lb_vv, yerr=lb_vv_std, fmt="*", color='blue', label="Data", markersize=6, elinewidth=1, capsize=1.5, barsabove=True)
        ax.errorbar(lb_sub, lb_vvc, fmt="o-", color='black', label="Calculated", markersize=2)
        ax.set_ylabel("1/V\u2080")
        ax.set_xlabel("1/[Substrate]")
        ax.set_title("Lineweaver-Burk Plot")

        plt.savefig(f'{name}.png', dpi=300)
        return fig
