function make_rundirectory(rundirname,bin_name)
% make a rundirectory with all the necessary files in the MODEL/runs directory
% this also compiles the MOREAD model
% the rundirectory will be under MODEL/runs
% Basically, MOREAD works in the following way:
% 1) Check that the FORTRAN CODE is like you want it
% 2) run make_rundirectory with the corret directory name and a descriptive
% binary name (bin_name)
% 3) cd to the rundirectory, and run run_MOREAD. It should recognise the
% binary as the only binary in the directory. Do not name other files with
% the suffix .bin

cd /Users/dal/MOREAD/MODEL/runs

eval(['!mkdir ' rundirname]);
eval(['!cp ../scripts/make_model_setup.m ' rundirname '/'])
eval(['!cp ../scripts/run_MOREAD.m ' rundirname '/'])
eval(['!cp ../scripts/format_fortran_d.m ' rundirname '/'])
eval(['!cp ../scripts/plot_disc_auto_exact.m ' rundirname '/'])
eval(['!cp ../scripts/disc_conv.m ' rundirname '/'])
eval(['!cp ../scripts/get_disc_nucrate_exact.m ' rundirname '/'])



% copy default sink file 
eval(['!cp ../inputs/SINKDIST_0.TXT ' rundirname '/'])

str = sprintf('cd %s',rundirname);
eval(str);

cd /Users/dal/MOREAD/MODEL/

% strip the .bin from the end if given, so no duplicate is made. 
% .bin is added automatically
if(strcmp(bin_name(end-3:end),'.bin'))
    bin_name=bin_name(end-3:end);
end

eval(['!/usr/local/bin/gfortran source/a_dvode_f90_m.f90 source/b_discrete_dvodef90.f90 source/discrete_helpers.f90 source/x_discrete.f90 -o binaries/' bin_name '.bin'])

eval(['!cp ./binaries/' bin_name '.bin ./runs/' rundirname '/'])

cd /Users/dal/MOREAD/MODEL/runs

str = sprintf('cd %s',rundirname);
eval(str);
