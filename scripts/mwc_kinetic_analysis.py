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

class MWC_Kinetic_Solver:
    """The equations and functions used to calculate the kinetic parameters for the Monod-Wyman-Changeux Model

    ...
    Usage
    -----
    MWC_Kinetic_Solver( v_t, v_r, k_t, k_r, l0, n )

    Attributes
    ----------
    V_T : int
        initial guess or solved velocity in tensed state value
    V_R : int
        initial guess or solved velocity in relaxed state value
    K_T : int
        initial guess or solved binding in tensed state value
    K_R : int
        initial guess or solved binding in relaxed state value
    L0 : int
        initial guess or solved for allosteric constant value
    n : int
        initial guess or solved number of allosteric sites
    """

    def __init__(self, V_T, V_R, K_T, K_R, L0, n):
        self.V_T = V_T
        self.V_R = V_R
        self.K_T = K_T
        self.K_R = K_R
        self.L0 = L0
        self.n = n
    
    def mwc_model(self, S):
        """MWC equation for a singular substrate index
        ...
        
        Attributes
        ----------
        S : list
            a list of substrate concentrations starting from most concentrated
        """

        V_T = self.V_T
        V_R = self.V_R
        K_T = self.K_T
        K_R = self.K_R
        L0 = self.L0
        n = self.n
        self.S = S    
        term_T = (1 + S / K_T) ** n
        term_R = (1 + S / K_R) ** n
        denom = L0 * term_T + term_R
        return (V_T * L0 * term_T + V_R * term_R) / denom

    def loss(self, params, S_data, v_data):
        """Compute the loss function (residual square sum) of the MWC fit

        Calculates the residiual sum of squares between the calculated and data velocity values.
        Adds a slight penalty for L0 and N parameters outside of typical biological constraints
        (guides the model to biological conditions unless fit is significantly better outside
        of these conditions) to prevent finding local minima.
        ...

        Attributes
        __________
        params = list
            list of the initial guesses for each MWC parameter
        S_data = list
            list of substrate dilution data
        v_data = lsit
            list of velocity values from enzyme kinetic assay data
        """

        V_T, V_R, K_T, K_R, L0, n = params
        ks = MWC_Kinetic_Solver(V_T, V_R, K_T, K_R, L0, n)
        v_model = ks.mwc_model(S_data)
        residual = np.sum((v_data - v_model) ** 2)
        if n < 1 or n > 8:
            penalty_n = 10 * (n - 8) ** 2  # quadratic penalty
        else:
            penalty_n = 0

        penalty_L0 = 0
        if L0 > 500:
            penalty_L0 = 100 * (L0 - 500) ** 2

        return residual + penalty_n + penalty_L0

    def square_sum(self, cal, dat):
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

    def sums(self, V_T, V_R, K_T, K_R, L0, n, vvalues, substrate):
        """Creates a total sum of the list of the square sum values for each substrate concentration
        ...
        Attributes
        ----------
        V_T : int
            initial guess or solved velocity in tensed state value
        V_R : int
            initial guess or solved velocity in relaxed state value
        K_T : int
            initial guess or solved binding in tensed state value
        K_R : int
            initial guess or solved binding in relaxed state value
        L0 : int
            initial guess or solved for allosteric constant value
        n : int
            initial guess or solved number of allosteric sites
        vvalues : list
            list of velocity values calculated from raw enzyme kinetic assay data
        substrate : list
            list of substrate dilutions
        """

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
        """Generates a symbolic representation of the full residual sum of squares equation using the MWC equation

        This function creates an extended version of the sum of squares equation above. The only values in the 
        full equation are the velocity values from data and the substrate concentration. Used later to find its
        global minimum to give V_T, V_R, K_T, K_R, L0, and n.
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
        """Calculates each partial derivative of the above extended sum of squares equation
        
        ...
        Attributes
        ----------
        equation : float
            equation generated from the full_equation function above
        """

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
        """Solves each partial derivative from above for each equal to zero

        This function sets each partial derivative of the full sum of squares equation
        to zero and finds a set of (h, vmax, km) that satisfies each partial equal to zero.
        This ensures that the solution found in the mimimum of the sum of squares equation, 
        determining the parameters that best fit the data. The initial MWC parameters guesses
        are used as starting guesses to speed up the calculation and ensure convergance.
        ...
        Attributes
        ----------
        df_dV_T : float
            partial derivative of the full equation, f, with respect to V_T
        df_dV_R : float
            partial derivative of the full equation, f, with respect to V_R
        df_dK_T : float
            partial derivative of the full equation, f, with respect to K_R
        df_dK_R : float
            partial derivative of the full equation, f, with respect to K_R
        df_dL0 : float
            partial derivative of the full equation, f, with respect to L0
        df_dn : float
            partial derivative of the full equation, f, with respect to N
        """

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
    """Import .csv file to compute velocity values and the linear range if applicable
    
    Reads CSV file and calculates the linear range using the find_first_and_last_small_diff
    function and the velocity values to determine the kinetic parameters using the 
    MWC_Kinetic_Solver class.
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
        return fig

    def cv_no_inset(self, cv_vval, bf_kp, cv_kp):
        """Graph both the best-fit and the cross-validation fit data
        ...
        Attributes
        ----------
        cv_vval : list
            list of the calculated velocity values from the Cross Validation method
        bf_kp : list
            list of the best-fit kinetic parameters
        cv_kp : list
            list of the cross validation kinetic parameters
        """

        substrate = self.substrate
        vvalues = self.vvalues
        vval_calc = self.vval_calc
        vv_std = self.vv_std
        self.cv_vval = cv_vval
        self.cv_kp = cv_kp
        self.bf_kp = bf_kp

        name = self.name
        fig, ax = plt.subplots(figsize=(4,3), dpi=250)
        plt.subplots_adjust(left=0.15, wspace=0.3, bottom=0.15)
        ax.errorbar(substrate, vvalues, yerr=vv_std, fmt="*", color='blue', label="Data", markersize=6, elinewidth=1, capsize=1.5, barsabove=True)
        ax.errorbar(substrate, vval_calc, fmt="o-", color='black', label="BF Calculated", markersize=2)
        ax.errorbar(substrate, cv_vval, fmt="o-", color='slateblue', label="CV Calculated", markersize=2)
        ax.set_ylabel("V\u2080")
        ax.set_xlabel("[Substrate]")
        ax.set_title("MWC Kinetic Plot")
        ax.set_xscale("log")
        hill_bf = '%.2f'%(bf_kp[5])
        hill_cv = '%.2f'%(cv_kp[5])
        ax.annotate(f'BF # of Allosteric Sites = {hill_bf}', 
            xy=(0.05, 0.97), # Adjust based on your bbox_to_anchor for the legend
            xycoords='axes fraction', 
            xytext=(0, 0), # Offset from xy in points
            textcoords='offset points',
            horizontalalignment='left', 
            verticalalignment='top',
            fontsize=7,
            fontstyle='italic')
        ax.annotate(f'CV # of Allosteric Sites = {hill_cv}', 
            xy=(0.05, 0.92), # Adjust based on your bbox_to_anchor for the legend
            xycoords='axes fraction', 
            xytext=(0, 0), # Offset from xy in points
            textcoords='offset points',
            horizontalalignment='left', 
            verticalalignment='top',
            fontsize=7,
            fontstyle='italic')
        fig.legend(loc='upper left', bbox_to_anchor=(0.01, 0.89), bbox_transform=plt.gca().transAxes, fontsize=8)
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
        ax.set_title("MWC Kinetic Plot")
        ax.set_xscale("log")
        fig.legend(loc='upper left', bbox_to_anchor=(0.01, 0.94), bbox_transform=plt.gca().transAxes, fontsize=8)
        ax.annotate(f'# Allosteric Sites = {hill}', 
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
            hill = "%.3f"%(kinetic_parameters_mut[i][5])
            ax.annotate(f'# Allosteric Sites = {hill}', 
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
        ax.set_title("MWC Kinetic Plot")
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

