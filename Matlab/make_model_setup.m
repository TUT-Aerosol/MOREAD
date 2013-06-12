function make_model_setup(in)


fid = fopen('MODEL_SETUP.TXT','w');

fprintf(fid,'IMAX=\n%i\n',in.imax);
fprintf(fid,'TEMP=\n%6.2f\n',in.temp);
fprintf(fid,'RH=\n%6.2f\n',in.rh);
fprintf(fid,'PRES=\n%6.2f\n',in.pres);
fprintf(fid,'NUCSIZE=\n%i\n',in.nucsize);
fprintf(fid,'NUCRATE=\n%s\n',format_fortran_d(in.nucrate));
fprintf(fid,'PULSE_LENGTH=\n%s\n',format_fortran_d(in.pulse_length));
fprintf(fid,'COND_ON=\n%i\n',in.cond_on);
fprintf(fid,'EVAP_ON=\n%i\n',in.evap_on);
fprintf(fid,'SINK_ON=\n%i\n',in.sink_on);
fprintf(fid,'COAG_ON=\n%i\n',in.coag_on);
fprintf(fid,'NUC_MECH=\n%i\n',in.nuc_mech);
fprintf(fid,'NUC_COEFF=\n%s\n',format_fortran_d(in.nuc_coeff));
fprintf(fid,'NUC_EXP=\n%s\n',format_fortran_d(in.nuc_exp));
fprintf(fid,'NUC_COEFF_ORG=\n%s\n',format_fortran_d(in.nuc_coeff_org));
fprintf(fid,'NUC_EXP_ORG=\n%s\n',format_fortran_d(in.nuc_exp_org));
fprintf(fid,'cvap_0=\n%s\n',format_fortran_d(in.cvap_0));

fclose(fid);