from calibration_functions import *

new_copied_res = np.array([  1.70396396e+02,   1.00000000e-01,   5.24324324e+04,
                             1.65535135e+04,   1.31594823e-01,   1.20969339e-01,
                             0.00000000e+00,   1.30891892e+02,   8.77371116e-01,
                             1.59224783e-01,   3.41341341e+01,   2.12939671e+00,
                             1.31173794e+00,   4.18799759e+00,   2.84775087e-01,
                             3.25698820e+00])
old_copied_res = np.array([ 1.67691692e+02,   4.63463504e-02,   9.90090090e+04,
                            1.55070270e+04,   6.76676625e-02,   1.46146135e-01,
                            7.77477470e+01,   2.74774758e+00,   2.21622388e+00,
                            5.12912874e-01,   5.98822777e+00,   2.79079058e-01,
                            8.02802797e+00])

in_data_loc = "../Data/calibration_data"
out_data = '../Data/final_calibration'
settings = pd.read_csv(in_data_loc+"/calibration_run_settings.tsv", index_col=0, sep="\t")
chem_df = pd.read_csv(in_data_loc+"/chem_data.tsv", sep="\t", index_col=[0,1], parse_dates=True)
gene_df = pd.read_csv(in_data_loc+"/gene_data.tsv", sep="\t", index_col=[0,1], parse_dates=True)
obs_data = chem_df.join(gene_df)

mod_types = [0, 2]
mod_res = [new_copied_res, old_copied_res]

these_settings = settings.copy()

mi, x = mod_types[0], mod_res[0]
mt = settings.columns[mi]
mb = settings.columns[mi+1]
to_optimize = list(settings.ix[settings.ix[:, mb], mt].index)
new_vals = pd.Series({i:j for i, j in zip(to_optimize, x)})
these_settings.loc[new_vals.index, mt] = new_vals
these_settings = these_settings.ix[settings.index, settings.columns]
arg_tuple = (these_settings.ix[:, mt], out_data+"/final/final_run.mat", obs_data.copy(), 'gene_objective')
subdf, out_f, obs_data, score_type = arg_tuple

if platform.system() == 'Linux':
  run_cmd = ""
else:
  run_cmd = "/Applications/MATLAB_R2016b.app/bin/"

model_loc = os.path.dirname(out_f)
source_dir = os.path.join(os.getcwd(), 'lake_model')
copyDirectory(source_dir, model_loc)
input_args_loc = os.path.join(model_loc, 'baseline.txt')
subdf.T.to_csv(input_args_loc, header=False, index=False, float_format='%g')
init_val_f = os.path.join(model_loc, "concs0.txt")
apply_conc_multiplier(subdf, init_val_f)
run_cmd = run_cmd +"matlab -nodisplay -nojvm -nosplash -nodesktop " 
run_cmd = run_cmd +"-r 'calibration_kaw; exit'"
output_loc = os.path.join(model_loc, 'outFile.txt')

with open(output_loc, 'w') as out_h:
    out_h.write(os.path.abspath(out_f))

p = sp.Popen(run_cmd, cwd=model_loc, shell=True, stderr=sp.PIPE, stdout=sp.PIPE)
stdout, stderr = p.communicate()
results_df = importratesandconcs_mod(out_f, 'full df')
return_scores = score_results(obs_data, results_df, "return_scores")

obs_df, data_df = obs_data.copy(), results_df.copy()
bool1 = data_df.index.isin(obs_df.index)
bool2 = obs_df.index.isin(data_df.index)
assert (bool2.sum() > 0) and (bool1.sum() > 0)
sub_data_df = data_df[bool1]
sub_obs_df = obs_df[bool2]
sub_data_df = combine_similar_procs(sub_data_df, sub_obs_df)
assert sub_obs_df.shape == sub_data_df.shape

plotdir= out_data+"/plots"
os.mkdir(plotdir)
for idx, col_ in enumerate(sub_obs_df.columns):
    obs_vec = sub_obs_df.ix[:, col_].dropna()
    mod_vec = sub_data_df.ix[:, col_].dropna()
    bool_11 = obs_vec.index.isin(mod_vec.index)
    bool_22 = mod_vec.index.isin(obs_vec.index)
    obs_vec_nn = obs_vec[bool_11]
    mod_vec_nn = mod_vec[bool_22]
    obs_vec_std, data_vec_std = standard_scale_df(obs_vec_nn), standard_scale_df(mod_vec_nn)
    obs_vec_std.name, data_vec_std.name = "observed", "modeled"
    col_pack = pd.concat((obs_vec_std, data_vec_std), 1)
    plt.clf()
    fig_, ax_ = plt.subplots(1, 1, figsize=(12, 14), dpi=300)
    col_pack.plot(ax=ax_, legend=True, rot=45, fontsize=14)
    plt.tight_layout()
    plt.savefig(plotdir+"/"+col_+".png", dpi=300)


#ts2 = these_settings.drop(["new_model_bool", "old_model_bool"], 1)
#score_df = pd.DataFrame(index=zip(*new_res_scores)[0], columns=ts2.columns)

