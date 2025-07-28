import math
import scipy as sc
import numpy as np
import pandas as pd
import argparse
from scipy.optimize import fsolve, least_squares, Bounds, minimize
from sympy import symbols, diff, solve, nsolve, checksol
import matplotlib
matplotlib.use('Agg')
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
    """Function to compare the slopes of kinetic data to determine linear range

    Calculates the percent differences between subsequent slope values to determine
    whether the difference is less than the threshold. If so, it uses the first 
    occurance as the 'first index' and the last occurance as the 'last index'.
    ...

    Attributes
    ----------
    data : list
        list of velocity values used to calculate percent differences from
    threshold : int
        the threshold with which to compare percent differences between subsequent slopes
        default is 10%
    ...
    """

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

class Hill_Kinetic_Solver:
    """The equations and functions used to calculate the kinetic parameters for the Hill Model

    ...
    Usage
    -----
    Hill_Kinetic_Solver( hill coefficient, v max value, km value)

    Attributes
    ----------
    h : int
        initial guess or solved Hill coefficient value
    vmax : int
        initial guess or solved Vmax value
    km : int
        initial guess or solved Km value
    """

    def __init__(self, h, vmax, km):
        self.h = h
        self.vmax = vmax
        self.km = km
        
    def hill_equation(self, num, substrate):
        """Hill equation for a singular substrate index
        ...
        
        Attributes
        ----------
        num : int
            index of substrate values to compute the Hill equation velocity for
        substrate : list
            a list of substrate concentrations starting from most concentrated
        """

        h = self.h 
        vmax = self.vmax
        km = self.km
        self.num = num
        self.substrate = substrate
        eq1 = (vmax*(substrate[num]**h))/((km**h)+(substrate[num]**h))
        return eq1
    
    def square_sum(self, dat, cal):
        """Calculates the square sum of a single calculated value and a data value
        ...
        Attributes
        ----------
        dat : int
            velocity value from enzyme kinetic assay data
        cal : int
            velocity value for a specific substrate calculated using the hill equation
        """

        x = cal
        y = dat
        eq2 = (x-y)**2
        return eq2

    def sums(self, h, vmax, km, vvalues, substrate):
        """Creates a total sum of the list of the square sum values for each substrate concentration
        ...
        Attributes
        ----------
        h : int
            the current guess for the Hill coefficient
        vmax : int
            the current guess for the Vmax value
        km : int
            the current guess for the Km value
        vvalues : list
            list of velocity values calculated from raw enzyme kinetic assay data
        substrate : list
            list of substrate dilutions
        """

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
        """Generates a symbolic representation of the full residual sum of squares equation using the Hill equation

        This function creates an extended version of the sum of squares equation above. The only values in the 
        full equation are the velocity values from data and the substrate concentration. Used later to find its
        global minimum to give h, vmax, and km.
        ...
        Attributes
        ----------
        substrate : list
            list of substrate dilutions
        vvalues : list
            list of velocity values from enzyme kinetic assay data
        """

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
        """Calculates each partial derivative of the above extended sum of squares equation
        
        ...
        Attributes
        ----------
        equation : float
            equation generated from the full_equation function above
        """

        h, vmax, km = symbols('h vmax km')
        f = equation
        df_dvmax = diff(f, vmax)
        df_dh = diff(f, h)
        df_dkm = diff(f, km)
        return df_dvmax, df_dh, df_dkm

    def minimize(self, df_dvmax, df_dh, df_dkm):
        """Solves each partial derivative from above for each equal to zero

        This function sets each partial derivative of the full sum of squares equation
        to zero and finds a set of (h, vmax, km) that satisfies each partial equal to zero.
        This ensures that the solution found in the mimimum of the sum of squares equation, 
        determining the parameters that best fit the data. The initial Hill coefficient, Vmax,
        and Km guesses are used as starting guesses to speed up the calculation and ensure 
        convergance.
        ...
        Attributes
        ----------
        df_dvmax : float
            partial derivative of the full equation, f, with respect to Vmax
        df_dh : float
            partial derivative of the full equation, f, with respect to Hill coefficient
        df_dkm : float
            partial derivative of the full equation, f, with respect to Km
        """

        h, vmax, km = symbols('h vmax km')
        eq1 = df_dvmax
        eq2 = df_dh
        eq3 = df_dkm
        sol = nsolve((eq1, eq2, eq3), (h, vmax, km), (self.h, self.vmax, self.km), prec=15, solver="bisect", verify=False)
        return [sol[0], sol[1], sol[2]]

class Import_Kinetic_Data:
    """Import .csv file to compute velocity values and the linear range if applicable
    
    Reads CSV file and calculates the linear range using the find_first_and_last_small_diff
    function and the velocity values to determine the kinetic parameters using the 
    Hill_Kinetic_Solver class.
    ...
    Attributes
    ----------
    fname : str
        file path to load CSV file from
    substrate : list
        list of substrate dilutions
    """

    def __init__(self, fname, substrate):
        self.fname = fname
        self.substrate = substrate
    
    def import_data(self, columns):
        """Reads CSV file and returns dataframe
        ...
        Attributes
        ----------
        columns : list
            list of first and last columns to read data from
        """

        self.columns = columns
        col_min = int(columns[0])
        col_max = int(columns[1])
        fname = self.fname
        df = pd.read_csv(fname, encoding='ISO-8859-1', usecols=range(col_min, col_max), nrows=31)
        df = df.apply(pd.to_numeric, errors="coerce")
        return df

    def gen_vvalues(self, df, time_min=15, time_max=40, steps=25, v_win=5):
        """Calculates the velocity values from the slopes of the raw kinetic data
        ...
        Attributes
        ----------
        df : dataframe
            pandas dataframe containing the raw kinetic assay data
        time_min : int
            the minimimum data point value to begin calculating slopes
            default is 15 (useful if data is taken every 1 min for 1 hour)
        time_max : int
            the maximum data point value to end calculating slopes
            defualt is 40
        steps : int
            the number of data points between time_min and time_max to calculate slopes from
            the value of steps must equally divise the difference between time_max and time_min
            default is 25
        v_win : int
            the window within which to calculate the slopes (a value of 5 means slopes 
            are calculated using (val5 - val1)/5)
            defualt is 5
        """

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
        """Function to estimate the linear range of the raw data

        Compares the slopes (calculated based on v_win) of subsequent 
        data points starting from the beginning of the data. Once the 
        percent difference is less that 2.5% (for the highest 4 substrate 
        concentrations) or 7.5% (for lower substrate concentrations), the
        function uses the latest first index from each substrate concentration
        and the earliest last index to determine the linear range that 
        encompases each substrate concentration.
        ...
        Attributes
        ----------
        df : dataframe
            pandas dataframe containing the raw enzyme kinetic assay data
        v_win : int
            window within which to calculate the slope
        """

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
        first_indices = [pair[0] for pair in lin_range_all if None not in pair]
        last_indices  = [pair[1] for pair in lin_range_all if None not in pair]
        avg_first = sum(first_indices) / len(first_indices)
        closest_first = min(first_indices, key=lambda x: abs(x - avg_first))
        avg_last = sum(last_indices) / len(last_indices)
        closest_last = min(last_indices, key=lambda x: abs(x - avg_last))
        return [closest_first, closest_last]


class get_inputs:
    """General class to generate some of the input values necessary to calculate the kinetic parameters
    """
    
    def find_nearest(self, array, value):
        """Finds the list value nearest to an input value
        ...
        Attributes
        ----------
        array : array
            any numpy array of interest to find the nearest value within
        value : int
            the value to find the most similar in the array
        """

        self.array = array
        self.value = value
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        return array[idx]

    def gen_substrate(self, sub):
        """Generate the list of substrate concentrations based on a serial dilution
        ...
        Attributes
        ----------
        sub : list
            list containing the number of substrate concentrations, the dilution factor, and the max substrate concentration
            in the form [# substrate, dilution factor, max substrate]
        """

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
        """Generate the x and y values corresponding to the linear hill equation
        
        Utilizes a normalized velocity value (normalized to the maximum velocity value)
        to take the logarithm.
        ...
        Attributes
        ----------
        vvalues : list
            list of velocity values calculated from the enzyme kinetic assay data
        substrate : list
            list of substrate dilutions
        """

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
        """Generates a linear regression line of the linearized Hill model
        ...
        Attributes
        ----------
        spx : list
            list of the logarithms of the substrate dilutions
        spy : list
            list of the logarithms of the normalized velocity values
        """

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
    """Class to graph the kinetic data generated with above functions
    ...
    Attributes
    ----------
    name : str
        name for the output file
    substrate : list
        list of substrate dilutions
    vvalues : list
        list of the velocity values from the enzyme kinetic assay data
    vval_calc : list
        list of the velocity values calculated based on the kinetic parameters matching the best fit data
    kinetic_parameters : list
        list of the best-fit kinetic_parameters: [hill coefficient, vmax, km]
    vv_std : list or int
        list of the standard deviations of the vvalues list if applicable, otherwise a value of 0 is used
    """

    def __init__(self, name, substrate, vvalues, vval_calc, kinetic_parameters, vv_std):
        self.name = name
        self.substrate = substrate
        self.vvalues = vvalues
        self.vval_calc = vval_calc
        self.kinetic_parameters = kinetic_parameters
        self.vv_std = vv_std

    def with_inset(self, spx, spy, linregx, poly1d_fn, bbox=(220, 150, 600, 500), ax_fs=5):
        """Graph of the Hill kinetic data using an inset plot of the linearized form of the Hill equation
        ...
        Attributes
        ----------
        spx : list
            list of the logarithms of the substrate dilutions
        spy : list
            list of the normalized logathrithms of the enzyme kinetic assay data
        linregx : list
            list of the linear regression of the spx and spy data
        poly1d_fn : obj
            object that defines the linear regression line
        bbox : tuple
            tuple describing the location of the inset graph
            default is (220, 150, 600, 500)
            *currently not used*
        ax_fs : int
            font size for the inset axis labels
            default is 5
        """

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
        ax.set_title("Hill Kinetic Plot")
        ax.set_xscale("log")
        inset_ax = inset_axes(ax, width="45%", height="35%", loc="upper left", bbox_to_anchor=(0.1, 0.04, 0.8, 0.9), bbox_transform=plt.gca().transAxes)
        inset_ax.set_xlabel("log[S]", fontsize=ax_fs)
        inset_ax.set_ylabel(r"$\log_ \frac{v}{(1-v)}$", fontsize=ax_fs, labelpad=-3)
        inset_ax.set_xscale("linear")
        inset_ax.plot(linregx, poly1d_fn(linregx), '--k')
        inset_ax.plot(spx, spy, ".", color='grey', markersize=4)
        inset_ax.tick_params(axis='both', labelsize=6)

        hill = '%.2f'%(kinetic_parameters[0])
        ax.annotate(f'Hill Coefficient = {hill}', 
            xy=(0.12, 0.97), # Adjust based on your bbox_to_anchor for the legend
            xycoords='axes fraction', 
            xytext=(0, 0), # Offset from xy in points
            textcoords='offset points',
            horizontalalignment='left', 
            verticalalignment='top',
            fontsize=7,
            fontstyle='italic')
        inset_ax.set_yticks(np.arange(-6, 4, 3))
        inset_ax.set_xticks(np.arange(-4, 5, 1))
        fig.savefig(f'{name}.png', dpi=300)
        fig.savefig(f'{name}.svg')
        return fig

    def no_inset(self):
        """Graph of the kinetic data with no inset plot
        """

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
        ax.set_title("Hill Kinetic Plot")
        ax.set_xscale("log")
        hill = '%.2f'%(kinetic_parameters[0])
        ax.annotate(f"Hill Coefficient = {hill}", xy=(5, 0), fontsize=7, fontstyle='italic')
        fig.savefig(f'{name}.png', dpi=300)
        fig.savefig(f'{name}.svg')
        return fig

    def rep_no_inset(self, vval_rep_avg, vval_calc_rep_avg, vval_rep_std, vval_calc_rep_std):
        """Graph replica data without inset plot
        ...
        Attributes
        ----------
        vval_rep_avg : list
            list of the averages of the replica velocity values from enzyme kinetic assay data
        vval_calc_rep_avg : list
            list of the averages of the replica calculated velocity values from each replica's best fit data
        vval_rep_std : list
            list of the standard deviations of replica velocity values
        vavl_calc_rep_std : list
            list of the standard deviations of replica calculated velocity values
        """

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
        ax.set_title("Hill Kinetic Plot")
        ax.set_xscale("log")
        fig.legend(loc='upper left', bbox_to_anchor=(0.01, 0.94), bbox_transform=plt.gca().transAxes, fontsize=7)
        ax.annotate(f'Hill Coefficient = {hill}', 
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
        """Generates a graph comparing Wild-Type and Mutant kinetic data
        ...
        Attributes
        ----------
        *note: vvalues and vval_calc from __init__ must be a list of lists containing both WT and Mutant data*
        kinetic_parameters_mut : list of lists
            a list of the list of best-fit kinetic parameters for WT and Mutant data
        """

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
            hill = "%.3f"%(kinetic_parameters_mut[i][0])
            ax.annotate(f'Hill Coefficient = {hill}', 
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
        ax.set_title("Hill Kinetic Plot")
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
        """Generates linear Lineweaver-Burk equation plot
        *note: currently not in use in GUI*
        """

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
