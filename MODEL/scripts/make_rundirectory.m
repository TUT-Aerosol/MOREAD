% make a rundirectory with all the necessary files in the MODEL/runs directory

rundirname = 'runs_slopepaper';

cd /Users/dal/MOREAD/MODEL/runs

eval(['!mkdir ' rundirname]);
eval(['!cp ../scripts/make_model_setup.m ' rundirname '/'])
eval(['!cp ../scripts/run_MOREAD.m ' rundirname '/'])
eval(['!cp ../scripts/format_fortran_d.m ' rundirname '/'])
eval(['!cp ../scripts/plot_disc_auto_exact.m ' rundirname '/'])
eval(['!cp ../scripts/disc_conv.m ' rundirname '/'])
eval(['!cp ../scripts/get_disc_nucrate_exact.m ' rundirname '/'])


eval(['!cp ../binaries/MOREAD_05.bin ' rundirname '/'])

% copy default sink file 
eval(['!cp ../inputs/SINKDIST_0.TXT ' rundirname '/'])

str = sprintf('cd %s',rundirname);
eval(str);