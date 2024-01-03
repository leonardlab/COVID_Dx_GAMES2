import os
from games.config.experimental_data import define_experimental_data
from games.utilities.saving import create_folder
from games.modules.sensitivity_analysis.sensitivity_analysis_chi2 import (
    all_param_sweeps_10pct
)
from games.plots.plots_sensitivity_analysis import tornado_plot_chi2

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
    sub_folder_name = "MODULE 4 - SENSITIVITY ANALYSIS"
    path = create_folder(folder_path, sub_folder_name)
    os.chdir(path)

    if "chi2" in settings["sensitivity_analysis"]:
        x, exp_data, exp_error = define_experimental_data(settings)
        pct_mse_low_list, pct_mse_high_list = all_param_sweeps_10pct(
            settings["parameters"], x, exp_data, exp_error, settings
        )
        tornado_plot_chi2(
            pct_mse_low_list, pct_mse_high_list,
            settings["free_parameter_labels"])