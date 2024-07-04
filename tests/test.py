import numpy as np
import pandas as pd
from SALib.sample import latin
from math import log10
import pickle


# [T7, RT, RNase, vRNA, cas13]
# [5.0, 2.5, 0.005, 1, 90]

#manual data conditions:
# [5.0, 2.5, 0.005, 0.0, 90]
# [5.0, 2.5, 0.005, 2.0, 90]
# [5.0, 2.5, 0.005, 20.0, 90]

#load manual data and process
# manual_data_path = "/Users/kdreyer/Library/CloudStorage/OneDrive-NorthwesternUniversity/KatieD_LL/Plant_Dx/Experimental_data/SARS_CoV2/Manual_data_8-1.xlsx"
# manual_data_df = pd.read_excel(manual_data_path, "Average_&_error")
# for col in manual_data_df:
#     label_str = manual_data_df.loc[0, col]
#     # print(label_str)
#     label_list = [float(i.strip(' []')) for i in label_str.split(',')]
#     manual_data_df.loc[0, col] = label_list
# # print(manual_data_df)
# manual_data_processed = manual_data_df.copy()[["Average", "Average.1", "Average.2"]]
# manual_error_processed = manual_data_df.copy()[["STDEV", "STDEV.1", "STDEV.2"]]
# print(manual_data_processed)
# print(manual_error_processed)

# out_path = "/Users/kdreyer/Documents/Github/COVID_Dx_GAMES2/src/games/config/"
# out_fname = "processed_manual_data_exp.pkl"
# out_fname_err = "processed_manual_data_err.pkl"

# manual_data_processed.to_pickle(out_path + out_fname)
# manual_error_processed.to_pickle(out_path + out_fname_err)

# df_manual = pd.read_pickle(out_path + out_fname)
# # print(df_manual)
# maxVal = 0.965130841737526
# for (columnName, columnData) in df_manual.items():
#     label = list(columnData.iloc[0])
#     data = list(columnData.iloc[1:])
#     print(label)
#     print(data)

# df_manual_err = pd.read_pickle(out_path + out_fname_err)
# # print(df_manual)
# maxVal = 0.965130841737526
# for (columnName, columnData) in df_manual_err.items():
#     label = list(columnData.iloc[0])
#     data = list(columnData.iloc[1:])
#     print(label)
#     print(data)


path_games = "/Users/kdreyer/Documents/Github/COVID_Dx_GAMES2/src/games/config/"
fname_csv = "initial_guesses_manual_rep1_init.csv"
fname_pkl = "initial_guesses_manual_rep1_init.pkl"
# with pd.option_context('display.precision', 10):
#     print(df)
# print(df.iloc[0]["k_cas13"])
# t = np.linspace(0, 240, (240*100)+1)

path = "/Users/kdreyer/Library/CloudStorage/OneDrive-NorthwesternUniversity/KatieD_LL/Plant_Dx/Plant-Dx_initial_guesses.xlsx"
# pd.set_option('display.float_format', '{:.9E}'.format)
# igs = pd.read_excel(path, sheet_name="Rep_1", dtype=np.float64)
# igs = igs.round(decimals=8)
# print(igs)
# igs_manual = igs.iloc[0:1]
# print(igs_manual)
# with pd.option_context('display.precision', 9):
#     print(igs)

# igs.to_pickle(path_games+fname_pkl)
# igs.to_csv(path_games+fname_csv)

# rep1_igs = pd.read_pickle(path_games+fname_pkl)
# with pd.option_context('display.precision', 11):
#     print(rep1_igs)
# print(rep2_igs)#.iloc[0]["k_cas13"])

# tspace = np.linspace(0, 120, (120 * 100) + 1)
# t = tspace[::300]
# # print(t, len(t))

# t2 = np.concatenate((np.linspace(0, 120, 41), np.linspace(0, 120, 41), np.linspace(0, 120, 41)))
# print(t2[82:])


param_names = ["k_cas13", "k_degv", "k_txn", 
               "k_FSS", "a_RHA", "b_RHA", "c_RHA", 
               "k_loc_deactivation", "k_scale_deactivation"]
param_vals = [[0.000224063, 0.000429181, 0.000524008], [136.0589787, 138.7508875, 153.6535126],
			  [1151.286829, 996.0075739, 970.3459417], [0.09603252, 0.050862836, 0.037781467], 
              [1.659857597, 1.328196638, 1.646849707], [12.483870364, 8.35144094, 10.63683507],
			  [55.712091373, 6.841518557, 34.5275701], [56.730991093, 58.57051396,59.69839426], 
              [6.906538862, 7.315648481, 7.561251886]]


																

ig_data = dict(zip(param_names, param_vals))

igs = pd.DataFrame(ig_data)

print(igs)
# igs.to_pickle(path_games+fname_pkl)

