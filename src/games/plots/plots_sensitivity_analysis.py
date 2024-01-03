import numpy as np
import matplotlib.pyplot as plt
from games.models.set_model import settings

plt.style.use(settings["context"] + "paper.mplstyle.py")

def tornado_plot_chi2(
        low_vals: list, high_vals: list, 
        param_labels: list
) -> None:

    """
    Creates a tornado plot for the sensitivity analysis

    Args:
        low_vals: a list of floats defining the percent changes 
            for decreasing each parameter by 10%

        high_vals: a list of floats defining the percent changes 
            for increasing each parameter by 10% 
        
        param_labels: a list of strings defining the parameter
            labels for the plot (Settings_COVID_Dx.py
            conditions_dictionary["real_param_labels_all"])

    
    Returns: none

    Figures:
        './tornado plot_' + model + '_' + data + tolerance + '.svg':
            tornado plot for the sensitivity analysis
    """
     
    num_params = len(param_labels)

    pos = np.arange(num_params) + .5 # bars centered on the y axis
    
    fig, (ax_left, ax_right) = plt.subplots(ncols=2)
    ax_left.set_title('Change in chi2 from mid to low', fontsize = 8)
    ax_right.set_title('Change in chi2 from mid to high', fontsize = 8)
    bars_left = ax_left.barh(pos, low_vals, align='center', facecolor='dimgrey')
    ax_right.set_yticks([])
    ax_left.bar_label(bars_left)
    ax_left.set_xlabel('% change in chi2')
    bars_right = ax_right.barh(pos, high_vals, align='center', facecolor='dimgrey')
    ax_left.set_yticks(pos)
    ax_right.bar_label(bars_right)
    ax_left.set_yticklabels(param_labels, ha='center', x=-.1)
    ax_right.set_xlabel('% change in chi2')
    # plt.show()
    plt.savefig('./tornado_plot_chi2.svg', dpi = 600)