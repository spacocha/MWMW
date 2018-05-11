% to run this file from the command line: `matlab -nodisplay -nojvm -nosplash -nodesktop -r 'calibration_kaw; exit'`

clear all

% read in the input data
concs0 = csvread(strcat(pwd, '/concs0.txt'));
% read in the default values (to be replaced dynamically)
param_vals = csvread(strcat(pwd, '/baseline.txt'));

% run the blasted thing
[time_slices, concs_history, rates_history] = lake( ...
	param_vals(1),	... oxygen_bubble_rate
	param_vals(2),	... nitrogen_source
	param_vals(3),	... nitrogen_ratio
	param_vals(4),	... carbon_source_from_5e4
	param_vals(5),	... oxygen_source
	param_vals(6),	... methane_source
        param_vals(7),	... t_max
	param_vals(8),	... sminus_precipitation
	param_vals(9),  ... fe_precipitation
	param_vals(10), ... carbon_precip
	param_vals(11),	... diffusion_constant
	param_vals(12), ... ma_op_o_fe_rate_const
	param_vals(13),	... ma_op_o_n_rate_const
	param_vals(14),	... ma_op_o_s_rate_const
	param_vals(15),	... ma_op_fe_n_rate_const
	param_vals(16), ... ma_op_s_n_rate_const
	param_vals(17), ... ma_op_ch4_o_rate_const
	param_vals(18), ... ma_op_ch4_n_rate_const
	param_vals(19),	... ma_op_ch4_s_rate_const
	param_vals(20),	... primary_ox_rate_const
	param_vals(21),	... c_lim_o
	param_vals(22),	... c_lim_n
	param_vals(23),	... c_lim_fe
	param_vals(24),	... c_lim_s
        concs0  ... Import from input file
);

% read in the file to pass out  
fileID = fopen(strcat(pwd,'/outFile.txt'),'r');
out_f = fscanf(fileID, '%s');
fclose(fileID);
% save the results
save(out_f, 'time_slices', 'rates_history', 'concs_history', '-v6');
