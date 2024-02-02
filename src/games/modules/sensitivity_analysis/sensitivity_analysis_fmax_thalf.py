import os
import numpy as np
from lmfit import Parameters, minimize
from copy import deepcopy
import json
from typing import Tuple
from alive_progress import alive_bar
from games.utilities.metrics import calc_r_sq, calc_percent_change
from games.models.set_model import model

def fitHill(y_exp: list) -> Tuple[float, float, float, float, float]:

    """
    Fits data to a hill function for a single set of
        conditions (component doses).

    Args: 
        y_exp: a list floats defining the normalized simulation 
            values for a single set of conditions

    Returns:
        f0: a float defining the initial value in the dataset

        fmax: a float defining the final value in the dataset

        km: a float defining the optimized parameter for t1/2

        n: a float defining the optimized parameter for the 
            Hill coefficient

        R_sq: a float defining the R squared value between the 
            data and the Hill fit
    """

    x = list(np.linspace(0, 240, 61)) #time (min)

    #Set v max to the final value of the time course
    fmax = y_exp[-1]
   
    #Set v0 to the intiial value of the time course
    f0 = y_exp[0]

    #Define a function to calculate the residual between the input simulation value (sim) and the Hill fit (model)
    def residual(p: list, x: list, y_exp: list) -> float:

        """
        Calculates the residual between the input simulation
            values and the Hill fit. Used in the minimization
            function as the cost function to be minimized
            between the simulation and the Hill fit.

        Args:
            p: a list of floats defining the parameters for
                the hill function

            x: a list of floats defining the time values for
                the simulation

            y_exp: a list of floats defining the simulation
                values

        Returns: 
            (y_exp - model): a float defining the residual 
            between the input simulation values and the Hill 
            fit
        """

        km = p['km'].value
        n = p['n'].value
        model = (((fmax - f0) * (x ** n)) / (km ** n + x ** n)) + f0
        return (y_exp - model)
    
    #Define parameters to be fit, along with initial guesses
    p = Parameters()
    p.add('km', value=10e-2, min=0, vary = True)
    p.add('n', value=2, min=0, max = 5, vary = True)
    
    #Perform the fit
    out = minimize(residual, p, args=(x, y_exp))
    bestFitParams = out.params.valuesdict()
    bestFitParamsList= bestFitParams.values()
    
    #Define fit parameters
    fit = []
    for value in bestFitParamsList:
        fit.append(value)
    km = fit[0]
    n = fit[1]
    
    #Simulate the Hill fit
    y_Hill = []
    for item in x:
        value = (((fmax - f0) * (item ** n)) / (km ** n + item ** n)) + f0
        y_Hill.append(value)
    
    #Calculate the R2 between the data and the Hill fit
    r_sq = calc_r_sq(y_exp, y_Hill)

    return f0, fmax, km, n, r_sq


def single_param_sweep_10pct(
        p: list, param_index: int, x: list[float],
        x_norm: list[float], settings: dict
) -> Tuple[float, float]:

    """
    Solves model ODEs for all conditions (component doses) for two 
        cases of changing parameter at param_index: increase by 10%
        and decrease by 10%

    Args:
        p: a list of floats defining the parameter values for all 
            potentially free parameters (Settings_COVID_Dx.py
            conditions_dictionary["p_all"])

        param_index: an integer defining the index of the parameter for the sweep

    Returns:
        mse_low: a float defining the mse resulting from the 10% decrease in
            the parameter

        mse_high: a float defining the mse resulting from the 10% increase in
            the parameter
    """

    p_vals = deepcopy(p)
    p_low = p_vals[param_index] * 0.9
    p_high = p_vals[param_index] * 1.1
    p_vals[param_index] = p_low
    model.parameters = deepcopy(p_vals)

    solutions_low, _ = model.solve_experiment([x])
    solutions_norm_condition_low, _ = model.solve_experiment([x_norm])

    if max(solutions_norm_condition_low) == 0:
        solutions_norm_low = [0]*len(solutions_low)
    else:
        solutions_norm_low = [i/max(solutions_norm_condition_low) for i in solutions_low]

    _, fmax_low, thalf_low, _, _ = fitHill(solutions_norm_low)

    p_vals[param_index] = p_high
    model.parameters = deepcopy(p_vals)

    solutions_high, _ = model.solve_experiment([x])
    solutions_norm_condition_high, _ = model.solve_experiment([x_norm])

    if max(solutions_norm_condition_high) == 0:
        solutions_norm_high = [0]*len(solutions_high)
    else:
        solutions_norm_high = [i/max(solutions_norm_condition_high) for i in solutions_high]

    _, fmax_high, thalf_high, _, _ = fitHill(solutions_norm_high)

   
    return fmax_low, thalf_low, fmax_high, thalf_high

def all_param_sweeps_fmax_thalf(
        p: list[float], parameter_labels: list[str], x: list[float],
        x_norm: list[float], settings: dict
 ) -> Tuple[list, list]:

    """
    Performs all parameter sweeps for increasing or decreasing 
        each parameter value by 10%

    Args:
        p: a list of floats defining the parameter values for all 
            potentially free parameters (Settings_COVID_Dx.py
            conditions_dictionary["p_all"])
    
    Returns:
        pct_mse_low_list: a list of floats defining the percent changes 
            for decreasing each parameter by 10%

        pct_mse_high_list: a list of floats defining the percent changes 
            for increasing each parameter by 10%
    """
    model.parameters = p
    solutions_mid, _ = model.solve_experiment([x])
    solutions_norm_condition_mid, _ = model.solve_experiment([x_norm])
    if max(solutions_norm_condition_mid) == 0:
        solutions_norm_mid = [0]*len(solutions_mid)
    else:
        solutions_norm_mid = [i/max(solutions_norm_condition_mid) for i in solutions_mid]

    _, fmax_mid, thalf_mid, _, _ = fitHill(solutions_norm_mid)

    # print('F_max mid: ', fmax_mid)
    # print('t_half mid: ', thalf_mid)

    
    pct_fmax_low_list = []
    pct_fmax_high_list = []
    pct_thalf_low_list = []
    pct_thalf_high_list = []
    with alive_bar(len(p)) as bar:
        for param_index in range(0, len(p)):
            (fmax_low, thalf_low,
             fmax_high, thalf_high) = single_param_sweep_10pct(
                 p, param_index, x, x_norm, settings
             )
            pct_fmax_low = calc_percent_change(fmax_mid, fmax_low)
            pct_fmax_low_list.append(pct_fmax_low)
            pct_fmax_high = calc_percent_change(fmax_mid, fmax_high)
            pct_fmax_high_list.append(pct_fmax_high)

            pct_thalf_low = calc_percent_change(thalf_mid, thalf_low)
            pct_thalf_low_list.append(pct_thalf_low)
            pct_thalf_high = calc_percent_change(thalf_mid, thalf_high)
            pct_thalf_high_list.append(pct_thalf_high)
            # print(
            #     parameter_labels[param_index],
            #     ": F_max low = ", fmax_low, "% Change F_max low = ",
            #     pct_fmax_low
            # )
            # print(
            #     parameter_labels[param_index],
            #     ": F_max high = ", fmax_high, "% Change F_max high = ",
            #     pct_fmax_high
            # )
            # print(
            #     parameter_labels[param_index],
            #     ": t_half low = ", thalf_low, "% Change t_half low = ",
            #     pct_thalf_low
            # )
            # print(
            #     parameter_labels[param_index],
            #     ": t_half high = ", thalf_high, "% Change t_half high = ",
            #     pct_thalf_high
            # )
            
            bar()

    return (
        pct_fmax_low_list, pct_fmax_high_list, 
        pct_thalf_low_list, pct_thalf_high_list
    )
