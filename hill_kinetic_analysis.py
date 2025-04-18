
#!/usr/bin/env python
# coding: utf-8

# In[1]:

import math
import scipy as sc
import numpy as np
import pandas as pd
import argparse
from scipy.optimize import fsolve, least_squares, Bounds, minimize
from sympy import symbols, diff, solve, nsolve, checksol
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# In[2]:

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


class Hill_Kinetic_Solver:
    """
    Kinetic_Solver( hill coefficient, v max value, km value)
    """

    def __init__(self, h, vmax, km):
        self.h = h
        self.vmax = vmax
        self.km = km
        
    def hill_equation(self, num, substrate):
        h = self.h 
        vmax = self.vmax
        km = self.km
        self.num = num
        self.substrate = substrate
        eq1 = (vmax*(substrate[num]**h))/((km**h)+(substrate[num]**h))
        return eq1
    
    def square_sum(self, cal, dat):
        x = cal
        y = dat
        eq2 = (x-y)**2
        return eq2

    def sums(self, h, vmax, km, vvalues, substrate):
        self.h = h
        self.vmax = vmax
        self.km = km
        self.vvalues = vvalues
        sums = []
        for i in range(len(vvalues)):
            ks = Hill_Kinetic_Solver(h, vmax, km)
            s = ks.square_sum(vvalues[i], ks.hill_equation(i, substrate))
            sums.append(s)
        value = sum(sums)
        return value

    def full_equation(self, substrate, vvalues):
        equation = 0
        self.substrate = substrate
        self.vvalues = vvalues 
        for i in range(len(substrate)):
            h, vmax, km = symbols('h vmax km')
            sub = substrate[i]
            vval = vvalues[i]
            eq1 = (vval - (vmax*(sub**h))/((km**h)+(sub**h)))**2
            equation += eq1
        return equation

    def partial_diff(self, equation):
        h, vmax, km = symbols('h vmax km')
        f = equation
        df_dvmax = diff(f, vmax)
        df_dh = diff(f, h)
        df_dkm = diff(f, km)
        return df_dvmax, df_dh, df_dkm

    def minimize(self, df_dvmax, df_dh, df_dkm):
        h, vmax, km = symbols('h vmax km')
        eq1 = df_dvmax
        eq2 = df_dh
        eq3 = df_dkm
        sol = nsolve((eq1, eq2, eq3), (h, vmax, km), (self.h, self.vmax, self.km), prec=15, solver="bisect", verify=False)
        return [sol[0], sol[1], sol[2]]



# In[1]:


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

    def gen_vvalues(self, df, time_min=15, time_max=40, steps=25):
        self.time_min = time_min
        self.time_max = time_max
        self.steps = steps
        substrate = self.substrate
        arr = df.to_numpy()
        vvalues = []
        for i in np.linspace(time_min, time_max, steps):
            vval_time = []
            for j in range(len(substrate)):
                i = int(i)
                k = int(i+10)
                v = abs(arr[i][j] - arr[k][j])/10
                vval_time.append(v)
            vvalues.append(vval_time)
        return vvalues


# In[ ]:

#def find_nearest(array, value):
#    array = np.asarray(array)
#    idx = (np.abs(array - value)).argmin()
#    return array[idx]

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
        for i in range(len(substrate)):
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

    def linreg(self, spx, spy, lin_min, lin_max):
        self.spx = spx
        self.spy = spy
        self.lin_min = lin_min
        self.lin_max = lin_max
        linregx = []
        linregy = []
        for i in range(6,17):
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
        ax.set_xlabel("[DHNP]")
        ax.set_title("Hill Kinetic Plot")
        ax.set_xscale("log")
        inset_ax = inset_axes(ax, width="45%", height="35%", loc=2, bbox_to_anchor=bbox)
        inset_ax.set_xlabel("log[S]", fontsize=ax_fs)
        inset_ax.set_ylabel(r"$\log_ \frac{v}{(1-v)}$", fontsize=ax_fs, labelpad=-3)
        inset_ax.set_xscale("linear")
        inset_ax.plot(linregx, poly1d_fn(linregx), '--k')
        inset_ax.plot(spx, spy, ".", color='grey', markersize=4)
        plt.yticks(fontsize=6)
        plt.xticks(fontsize=6)

        hill = '%.2f'%(kinetic_parameters[0])
        ax.annotate(f"Hill Coefficient = {hill}", xy=(5, 0), fontsize=7, fontstyle='italic')
        plt.yticks(np.arange(-6, 4, 3))
        plt.xticks(np.arange(-4, 5, 1))
        plt.savefig(f'{name}.png')
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
        ax.set_xlabel("[DHNP]")
        ax.set_title("Hill Kinetic Plot")
        ax.set_xscale("log")
        hill = '%.2f'%(kinetic_parameters[0])
        ax.annotate(f"Hill Coefficient = {hill}", xy=(5, 0), fontsize=7, fontstyle='italic')
        plt.savefig(f'{name}.png')
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
        ax.set_xlabel("1/[DHNP]")
        ax.set_title("Lineweaver-Burk Plot")

        plt.savefig(f'{name}.png')
        return fig
