function[in_updated] = make_model_setup(in)
% sets up the model:
% MODEL_SETUP.TXT
% 
in_updated = in;

fprintf('Writing MODEL_SETUP...')
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
fprintf(fid,'Qvap_0=\n%s\n',format_fortran_d(in.qvap_0));
fprintf(fid,'CONST_CVAP=\n%i\n',in.const_cvap);


fclose(fid);
fprintf('.done.\n')


% creating the the sink file and updating the input structure with the sink
% distribution
csfile = textread(in.sinkfilename); 
cs_size = csfile(1,1);
CS_Dp = csfile(2,1:cs_size).*2;
CS_N  = csfile(3,1:cs_size);

CS_0 = CS_general(CS_Dp,CS_N, in.temp,1.0)
CSscale = in.condsink_value./CS_0; % the scaling factor
CS_N_run = CSscale.*CS_N;

% write the new CS file
if exist('SINKDIST.TXT')
    delete('SINKDIST.TXT')
end
fprintf('Copying and scaling sink file %s to SINKDIST.TXT...',in.sinkfilename)
fid = fopen('SINKDIST.TXT','w');
fprintf(fid,'%i\n',cs_size); % the distribution length
fprintf(fid,'%s\n',format_fortran_d(CS_Dp./2)); % scale back to diameter
fprintf(fid,'%s\n',format_fortran_d(CS_N_run));
fclose(fid);

in_updated.sinkdist = [CS_Dp(:)'; CS_N_run(:)'];


fprintf('.done.\n')

in_updated