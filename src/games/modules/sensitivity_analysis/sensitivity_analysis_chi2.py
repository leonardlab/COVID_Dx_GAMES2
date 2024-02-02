import os
import numpy as np
from copy import deepcopy
from typing import Tuple
from alive_progress import alive_bar
from games.modules.solve_single import solve_single_parameter_set
from games.utilities.metrics import calc_percent_change
from games.models.set_model import model

def single_param_sweep_10pct(
        p: list, param_index: int, x: list[float],
        exp_data: list[float], exp_error: list[float],
        settings: dict
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
    _, _, mse_low, _ = solve_single_parameter_set(
        x, exp_data, exp_error, settings["weight_by_error"]
    )

    p_vals[param_index] = p_high
    model.parameters = deepcopy(p_vals)
    _, _, mse_high, _ = solve_single_parameter_set(
        x, exp_data, exp_error, settings["weight_by_error"]
    )   
    return mse_low, mse_high


def all_param_sweeps_mse(
        p: list[float], parameter_labels: list[str], x: list[float],
        exp_data: list[float], exp_error: list[float],
        settings: dict
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
    _, _, mse_mid, _ = solve_single_parameter_set(
        x, exp_data, exp_error, settings["weight_by_error"]
    )
    # print('mse opt: ', mse_mid)
    
    pct_mse_low_list = []
    pct_mse_high_list = []
    with alive_bar(len(p)) as bar:
        for param_index in range(0, len(p)):
            mse_low, mse_high = single_param_sweep_10pct(
                p, param_index, x, exp_data, exp_error, settings
            )
            pct_mse_low = calc_percent_change(mse_mid, mse_low)
            pct_mse_low_list.append(pct_mse_low)
            pct_mse_high = calc_percent_change(mse_mid, mse_high)
            pct_mse_high_list.append(pct_mse_high)
            # print(
            #     parameter_labels[param_index],
            #     ": MSE low = ", mse_low, "% Change MSE low = ",
            #     pct_mse_low
            # )
            # print(
            #     parameter_labels[param_index],
            #     ": MSE high = ", mse_high, "% Change MSE high = ",
            #     pct_mse_high
            # )
            bar()

    return pct_mse_low_list, pct_mse_high_list


