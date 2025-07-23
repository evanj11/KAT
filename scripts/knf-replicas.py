import knf_kinetic_analysis
import numpy as np
import pandas as pd
import math
import ast
import os
from numpy.random import default_rng
from sklearn.utils import resample
from sklearn.model_selection import KFold, ShuffleSplit   
from scipy.optimize import minimize, differential_evolution
from knf_kinetic_analysis import get_parser, Import_Kinetic_Data, KNF_Kinetic_Solver, get_inputs, graph_kinetic_data


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

def bayesian_bootstrap(sites, s, v, best_params, n_iter=1000):
    rng = default_rng(seed=42)
    s, v = np.array(s), np.array(v)
    n = len(s)
    boot_params = []
    bounds = [
        (0, 10000),       # V_T
        (0, 10000),      # V_R
        (0, 10000),   # K_T
        (-50, 50),    # K_R
    ]
    val = KNF_Kinetic_Solver(sites, best_params[0], best_params[1], best_params[2], best_params[3])
    for _ in range(n_iter):
        weights = rng.dirichlet(np.ones(n)*4)
        def weighted_rss(params):
            v_pred = val.knf_model(s)
            return np.sum(weights * (v - v_pred)**2)
        for _ in range(5):
            guess = best_params * rng.uniform(0.5, 1.5, size=len(best_params))
            result = minimize(weighted_rss, guess, bounds=bounds)
        if result.success:
            boot_params.append(result.x)

    return np.array(boot_params)

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
    conf_all = []
    for i in range(len(pathds)):
        data = Import_Kinetic_Data(pathds[i], substrate)
        df_data = data.import_data(columns)
        if len(pathds) == len(pathbs):
            blank = Import_Kinetic_Data(pathbs[i], substrate)
            df_blank = blank.import_data(columns)
            df = pd.DataFrame(df_data.values - df_blank.values)
            print(df)
        elif 1 == len(pathbs):
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
        
        with open(os.path.join(work_dir, 'sites.txt'), 'r') as file:
            lines = [line.strip() for line in file.readlines()]
            print(lines)
            if len(lines[0].split(",")) > 1:
                site_all = lines[0].split(",")
                print(site_all)
                sites = site_all[i]
            else:
                sites = lines[0]


        print(vv_std)

        vmax_guess = vvalues[0]
        kd_guess = substrate[np.abs(vvalues[0] - 0.5 * vmax_guess).argmin()]
        kbasal_guess = np.min(vvalues)
        val = KNF_Kinetic_Solver(sites, vmax_guess, kd_guess, kbasal_guess, 2)
        
        p0 = [vmax_guess, kd_guess, kbasal_guess, 2.0]

        bounds = [(0, 100000), (0, 10000), (0, 10000), (-50, 50)]

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

       # Optimized parameters
        kinetic_parameters = refined.x
        print("Best fit parameters:", kinetic_parameters)
        print("Minimized residual sum of squares:", refined.fun)       

        X = np.array(substrate)
        y = np.array(vvalues)

        kf = KFold(n_splits=5, shuffle=True, random_state=42)
        cv_params = []
        errors = []

        for train_idx, test_idx in kf.split(X):
            X_train, y_train = X[train_idx], y[train_idx]
            kf_res = differential_evolution(
                val.loss, 
                args=(X_train, y_train),
                bounds=bounds,
                strategy='best1bin',
                maxiter=1000,
                popsize=20,
                tol=1e-7
                )
            kf_ref = minimize(val.loss, kf_res.x, args=(X_train, y_train), bounds=bounds, method='L-BFGS-B', options={
                'maxiter': 10000,
                'ftol': 1e-12,    # tighter convergence on function value
                'gtol': 1e-10,    # tighter gradient convergence
                'eps': 1e-8       # step size
            })
            if kf_res.success:
                cv_params.append(kf_res.x)
                errors.append(kf_res.fun)
        cv_mean_params = np.mean(cv_params, axis=0)
        errors_mean = np.mean(errors, axis=0)
        print(f"Params = {cv_mean_params}")
        print(f"RSS = {errors_mean}")

        multi_start_result = None
        multi_start_loss = np.inf
    
        for _ in range(100):
            guess = np.random.uniform([0, 0, 0, -50], [1000, 1000, 1000, 50])
            res = minimize(val.loss, guess, args=(substrate, vvalues), bounds=bounds, method='L-BFGS-B', options={
                'maxiter': 10000,
                'ftol': 1e-12,
                'gtol': 1e-10,
                'eps': 1e-8})
            if res.fun < multi_start_loss:
                multi_start_loss = res.fun
                multi_start_result = res.x
        print(multi_start_result)

        if len(substrate) < 30:
            # Begin bootstrapping data
            n_bootstrap = 1000
            rng = default_rng(seed=42)
            boot_params = []

            boot_params = bayesian_bootstrap(sites, substrate, vvalues, cv_mean_params, n_iter=1000)
       
            boot_params = np.array(boot_params)
            param_names = ['Vmax', 'Kd', 'K_basal', 'Gamma']
            ci_lower = np.percentile(boot_params, 0.5, axis=0)
            ci_upper = np.percentile(boot_params, 99.5, axis=0)
            param_names = ["V_T", "V_R", "K_T", "K_R", "L0", "n"]
            for names, low, high in zip(param_names, ci_lower, ci_upper):
                print(f"{names}: 99% CI = [{low:.3f}, {high:.3f}]")


        #AIC and BIC checks
        n = len(substrate)
        aic_bf = n * math.log(refined.fun / n) + 6 * 2
        aic_kf = n * math.log(errors_mean / n) + 6 * 2
        aic_diff = abs(aic_bf - aic_kf)

        bic_bf = n * math.log(refined.fun / n) + 6 * math.log(n)
        bic_kf = n * math.log(errors_mean / n) + 6 * math.log(n)
        bic_diff = bic_bf - bic_kf
        print(f'BIC results = {bic_diff}')
        if bic_diff < 0:
            print('Best Fit has lower RSS')
        else:
            print('KFold has lower RSS')

        if len(substrate) < 30:
            within_ci = []
            for names, best, low, high in zip(param_names, kinetic_parameters, ci_lower, ci_upper):
                if low <= best <= high:
                    within_ci.append('True')
                else:
                    within_ci.append('False')

        confidence = []
        for val in range(4):
            val_diff = (multi_start_result[val] - refined.x[val]) / refined.x[val]
            if val_diff < 0.05:
                confidence.append('True')
            else:
                confidence.append('False')



        val = KNF_Kinetic_Solver(sites, kinetic_parameters[0], kinetic_parameters[1], kinetic_parameters[2], kinetic_parameters[3])
        vval_calc = val.knf_model(substrate)

        with open(os.path.join(work_dir, 'name_data.txt'), 'r') as file:
            file_name = [line.strip() for line in file.readlines()]
            print(file_name)

        with open(os.path.join(work_dir, f'{file_name[0]}_{i}.txt'), 'w') as file:
            file.write('Complex model used')
            file.write("\n")
            file.write('Koshland-Nemethy-Filmer')
            file.write("\n")
            file.write('Vmax = ' + str(kinetic_parameters[0]))
            file.write("\n")
            file.write('Kd = ' + str(kinetic_parameters[1]))
            file.write("\n")
            file.write('K_Basal = ' + str(kinetic_parameters[2]))
            file.write("\n")
            file.write('Gamma = ' + str(kinetic_parameters[3]))
            file.write("\n")
            file.write('Minimized residual sum of squares: ' + str(refined.fun))
            file.write("\n\n")
            file.write('Cross-Validation (KFold, k=5) Mean Parameters:\n')
            file.write('Vmax = ' + str(cv_mean_params[0]))
            file.write("\n")
            file.write('Kd = ' + str(cv_mean_params[1]))
            file.write("\n")
            file.write('K_Basal = ' + str(cv_mean_params[2]))
            file.write("\n")
            file.write('Gamma = ' + str(cv_mean_params[3]))
            file.write("\n")
            file.write('Minimized residual sum of squares: ' + str(errors_mean))
            file.write("\n\n")
            if len(substrate) < 30:
                for names, low, high in zip(param_names, ci_lower, ci_upper):
                    file.write(f"{names}: 99% CI = [{low:.3f}, {high:.3f}]\n")
                if 'False' in within_ci:
                    file.write('\n')
                    file.write('There is a high possibility for overfitting, use Cross-Validation values')
                else:
                    file.write('\n')
                    file.write('Each parameter is within the 99% confidence interval range')
        if len(substrate) < 30:
            if 'False' in within_ci:
                kinetic_parameters = cv_mean_params
            else:
                kinetic_parameters = refined.x

        with open(os.path.join(work_dir, f'complex_data_{i}.txt'), 'w') as file:
            file.write('Complex model used')
            file.write("\n")
            file.write('Koshland-Nemethy-Filmer')
            file.write("\n")
            file.write(str(kinetic_parameters[0]))
            file.write("\n")
            file.write(str(kinetic_parameters[1]))
            file.write("\n")
            file.write(str(kinetic_parameters[2]))
            file.write("\n")
            file.write(str(kinetic_parameters[3]))
            file.write("\n")
            if len(substrate) < 30:
                if 'False' in within_ci:
                    file.write('\n')
                    file.write('Poor confidence\n')
                    file.write('Using Cross-Validation values')

        vvalues_rep.append(vvalues)
        vval_calc_rep.append(vval_calc)
        kinetic_parameters_rep.append(kinetic_parameters)
        if len(substrate) < 30:
            conf_all.extend(within_ci)
        plot = graph_kinetic_data(os.path.join(work_dir, f"{file_name[0]}_{i}"), substrate, vvalues, vval_calc, kinetic_parameters, 0)
        plot.no_inset()
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
    for m in range(4):
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

    with open(os.path.join(work_dir, f'{file_name[0]}.txt'), 'w') as file:
        file.write('Complex model used')
        file.write('\n')
        file.write('Koshland-Nemethy-Filmer')
        file.write("\n")
        vmax_avg = '%.3f'%(kinetic_parameters_avg[0])
        vmax_std = '%.3f'%(kinetic_parameters_std[0])
        file.write('Vmax = ' + str(vmax_avg))
        file.write(" \u00B1 ")
        file.write(str(vmax_std))
        file.write("\n")
        kd_avg = '%.3f'%(kinetic_parameters_avg[1])
        kd_std = '%.3f'%(kinetic_parameters_std[1])
        file.write('Kd = ' + str(kd_avg))
        file.write(" \u00B1 ")
        file.write(str(kd_std))
        file.write("\n")
        kb_avg = '%.3f'%(kinetic_parameters_avg[2])
        kb_std = '%.3f'%(kinetic_parameters_std[2])
        file.write('K_Basal = ' + str(kb_avg))
        file.write(" \u00B1 ")
        file.write(str(kb_std))
        file.write("\n")
        y_avg = '%.3f'%(kinetic_parameters_avg[3])
        y_std = '%.3f'%(kinetic_parameters_std[3])
        file.write('Gamma = ' + str(y_avg))
        file.write(" \u00B1 ")
        file.write(str(y_std))
        file.write("\n")
        
    with open(os.path.join(work_dir, f'{file_name[0]}_complex_data.txt'), 'w') as file:
        file.write('Complex model used')
        file.write('\n')
        file.write('Koshland-Nemethy-Filmer')
        file.write("\n")
        vmax_avg = '%.3f'%(kinetic_parameters_avg[0])
        vmax_std = '%.3f'%(kinetic_parameters_std[0])
        file.write(str(vmax_avg))
        file.write(" \u00B1 ")
        file.write(str(vmax_std))
        file.write("\n")
        kd_avg = '%.3f'%(kinetic_parameters_avg[1])
        kd_std = '%.3f'%(kinetic_parameters_std[1])
        file.write(str(kd_avg))
        file.write(" \u00B1 ")
        file.write(str(kd_std))
        file.write("\n")
        kb_avg = '%.3f'%(kinetic_parameters_avg[2])
        kb_std = '%.3f'%(kinetic_parameters_std[2])
        file.write(str(kb_avg))
        file.write(" \u00B1 ")
        file.write(str(kb_std))
        file.write("\n")
        y_avg = '%.3f'%(kinetic_parameters_avg[3])
        y_std = '%.3f'%(kinetic_parameters_std[3])
        file.write(str(y_avg))
        file.write(" \u00B1 ")
        file.write(str(y_std))
        file.write("\n")
        if len(substrate) < 30:
            if 'False' in conf_all:
                file.write('\n')
                file.write('Poor confidence\n')
                file.write('Using Cross-Validation values')
       
    y_tot = y_avg + '\u00B1' + y_std
    plot_rep = graph_kinetic_data(os.path.join(work_dir, file_name[0]), substrate, vvalues_rep, vval_calc_rep, y_tot, 0)
    with open(os.path.join(work_dir, "mutant.txt"), 'r') as file:
        mut = [line.strip() for line in file.readlines()]

    if 'mutant' in mut:
        with open(os.path.join(work_dir, f'{file_name[0]}_mutant_data.txt'), 'w') as file:
            for mut in range(len(kinetic_parameters_rep)):
                file.write(str("%.3f"%(kinetic_parameters_rep[mut][0])))
                file.write(",")
                file.write(str("%.3f"%(kinetic_parameters_rep[mut][1])))
                file.write(",")
                file.write(str("%.3f"%(kinetic_parameters_rep[mut][2])))
                file.write(",")
                file.write(str("%.3f"%(kinetic_parameters_rep[mut][3])))
                file.write("\n")
        plot_rep.mut_rep(kinetic_parameters_rep)
    else:
        plot_rep.rep_no_inset(vval_rep_avg, vval_calc_rep_avg, vval_rep_std, vval_calc_rep_std)

if __name__ == "__main__":
    main()

