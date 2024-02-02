import os
from copy import deepcopy
from games.config.experimental_data import define_experimental_data
from games.utilities.saving import create_folder
from games.modules.sensitivity_analysis.sensitivity_analysis_chi2 import (
    all_param_sweeps_mse
)
from games.modules.sensitivity_analysis.sensitivity_analysis_fmax_thalf import (
    all_param_sweeps_fmax_thalf
)
from games.plots.plots_sensitivity_analysis import tornado_plot

def run_sensitivity_analysis(
        settings:dict, folder_path: str
) -> None:
    """Runs sensitivity analysis specified in settings

    Parameters
    ----------
    settings
        a dictionary of run settings

    folder_path
        a string defining the path to the main results folder

    Returns
    -------
    None
    
    """
    sub_folder_name = "MODULE 4 - SENSITIVITY ANALYSIS " + settings["dataID"][:4].upper()
    path = create_folder(folder_path, sub_folder_name)
    os.chdir(path)

    if "chi2" in settings["sensitivity_analysis"]:
        if settings["mechanismID"] == "D":
            p = deepcopy(settings["parameters "+settings["dataID"][:4]])
            parameter_labels = settings["free_parameter_labels"]
        elif settings["mechanismID"] == "D separate free params and fixed params":
            p = deepcopy(settings["parameters "+settings["dataID"][:4]] + settings["fixed_parameters"])
            parameter_labels = settings["free_parameter_labels"] + settings["fixed_parameter_labels"]

        x, exp_data, exp_error = define_experimental_data(settings)
        pct_mse_low_list, pct_mse_high_list = all_param_sweeps_mse(
            p, parameter_labels, x, exp_data, exp_error, settings
        )
        tornado_plot(
            pct_mse_low_list, pct_mse_high_list,
            "MSE", parameter_labels
        )
    if "f_max and t_half" in settings["sensitivity_analysis"]:
        if settings["mechanismID"] == "D":
            p = deepcopy(settings["parameters "+settings["dataID"][:4]])
            parameter_labels = settings["free_parameter_labels"]
        elif settings["mechanismID"] == "D separate free params and fixed params":
            p = deepcopy(settings["parameters "+settings["dataID"][:4]] + settings["fixed_parameters"])
            parameter_labels = settings["free_parameter_labels"] + settings["fixed_parameter_labels"]
        
        x_mid = [5.0, 2.5, 0.005, 1, 90]
        x_norm = [1.0, 2.5, 0.005, 10, 90]

        (pct_fmax_low_list, pct_fmax_high_list, 
         pct_thalf_low_list, pct_thalf_high_list
        ) = all_param_sweeps_fmax_thalf(p, parameter_labels, x_mid, x_norm, settings)

        tornado_plot(
            pct_fmax_low_list, pct_fmax_high_list,
            "F_max", parameter_labels
        )
        tornado_plot(
            pct_thalf_low_list, pct_thalf_high_list,
            "t_half", parameter_labels
        )

        