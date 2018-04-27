import pandas as pd
import cPickle as pickle

f_res_loc = "../Data/calibration_results/final_result.p"
with open(f_res_loc, "rb") as p_fh:
   dump_d = pickle.load(p_fh)

param_trace, dev_trace = dump_d
writer = pd.ExcelWriter('../Data/calibration_results/model_comparison.xlsx', engine='xlsxwriter')

sheet1 = (param_trace['old_model'].ix[:, 14] - param_trace['new_model'].ix[:, 17]).to_frame(name="Param_Val_Diff")
sheet1["old_model_start"] = param_trace['old_model'].ix[:, 0]
sheet1["old_model_end"] = param_trace['old_model'].ix[:, 14]
sheet1["new_model_start"] = param_trace['new_model'].ix[:, 0]
sheet1["new_model_end"] = param_trace['new_model'].ix[:, 17]
sheet1.to_excel(writer, sheet_name='Param Value Results')

pt_df = param_trace["new_model"]
mean_param_val_vec = pt_df.values.mean(1)
dpdx = (abs(pt_df.values[:, :-1] - pt_df.values[:, 1:]) / mean_param_val_vec[:, None])
dpdx_df = pd.DataFrame(data=dpdx, index=pt_df.index)
dpdx_df.to_excel(writer, sheet_name="new_mod_change_per_turn")

ave_param = param_trace["new_model"].mean(1)[:-1].values
param_std = dev_trace["new_model"].fillna(0).values
param_cvs = param_std/ave_param[:, None]
new_cvs = pd.DataFrame(index=dev_trace["new_model"].index, data=param_cvs)
new_cvs.to_excel(writer, sheet_name="new_mod_sampling_cv")

pt_df = param_trace["old_model"]
mean_param_val_vec = pt_df.values.mean(1)
dpdx = (abs(pt_df.values[:, :-1] - pt_df.values[:, 1:]) / mean_param_val_vec[:, None])
dpdx_df = pd.DataFrame(data=dpdx, index=pt_df.index)
dpdx_df.to_excel(writer, sheet_name="old_mod_change_per_turn")

ave_param = param_trace["old_model"].mean(1)[:-1].values
param_std = dev_trace["old_model"].fillna(0).values
param_cvs = param_std/ave_param[:, None]
old_cvs = pd.DataFrame(index=dev_trace["old_model"].index, data=param_cvs)
old_cvs.to_excel(writer, sheet_name="old_mod_sampling_cv")
writer.save()


