function[out,in] = run_MOREAD(varargin)
% this initializes the MOREAD model and sets the parameters, and then runs
% the model
% if no arguments are given, the default values will be used

binaryNameGiven =0;


% defining the deafult values
 in.imax = 1000;
 in.temp = 273.15;
 in.rh = 0.75;
 in.pres = 1.0;
 in.nucsize = 80;
 in.nucrate = 1000.*1e6;   % only if nuc_mech = 3 (constant)
 in.pulse_length = 1800;    
 in.cond_on = 1;
 in.evap_on = 0;
 in.sink_on = 1;
 in.coag_on = 1;
 in.nuc_mech = 1; % 1 = act 2 = kin 3 = const, 4 = free exp else = 0
 in.nuc_coeff = 1e-14; % needs to correspond to the nuc_mech
 in.nuc_exp = 3;  % only if nuc-mech = 4
 in.nuc_coeff_org = 1e-14; % not used
 in.nuc_exp_org = 2; % not used
 in.cvap_0 = 1e7.*1e6;  % starting cvap
 in.qvap_0 = 1e5.*1e6;  % starting vapour source
 in.const_cvap = 1; 
 in.sinkfilename = 'SINKDIST_0.TXT'; % this is the CS distribution shape file
 in.runfilename  = 'TEST_MOREAD';
 in.rundate = datestr(now);
 in.model_version = 0.5;
 in.binary_name = 'MOREAD_05.bin';
 in.condsink_value = 1e-3;  % the CS distribution will be scaled so that the CS matches this value

param_names = fieldnames(in);
total_params = length(param_names);


% If user has set parameters in the input, replace default values with
% user-defined values.

% Check that the input is even (e.g. initialize('gas_source',1e5)):
if(rem((nargin),2) ~= 0)
    error('set_initials: The number of arguments must be even.');
end

set_parameters=(nargin)/2;    % Get the number of arguments set in the input

% Check for duplicate definitions:
for i=1:2:set_parameters*2
    summary = 0;
    for j=1:2:set_parameters*2
        summary=summary+strcmp(varargin{i},varargin{j});
        if(summary > 1)
            error('set_initials: Parameter %s is defined twice.', varargin{i});
        end
    end
end



for i=1:2:set_parameters*2  % Go through all parameter names in the input
      parameter_name=varargin{i};
      if(strcmp(parameter_name,'binary_name'))
            % if the parameter is the binary name, set a flag that it has
            % been given
            binaryNameGiven = 1;
      end
      
      for j=1:total_params % Go through all param fields
          % Compare the input name with the names of param fields:
          if(strcmp(parameter_name,param_names{j}))
              % If the names match, replace the default param value with 
              % input value:
              in.(parameter_name)=varargin{i+1};
              break;
          elseif(j==total_params)
              % If we are in the end of loop and no matching name is found,
              % such a field doesn't exist.
              warning('set_initials: Invalid argument: %s', parameter_name);              
          end
      end
end

% if no binary name is given in the arguments,
% find the binary
if ~binaryNameGiven
    bn = dir('*.bin');
    if length(bn)>1,
        error('Multiple binary files found in rundirectory. Please give the binary file as an input value. Example: run_MOREAD(''binary_name'',''MOREAD.bin'')')
    elseif isempty(bn),
        error('No binary found. Compile or copy a binary file and place in the current directory.')
    elseif length(bn)==1,
        fprintf('Using only binary in directory: %s/n',bn.name);
        in.binary_name = bn.name;
    else
        error('Something really strange happened. Get a new computer.')
    end
end
 
 setenv('DYLD_LIBRARY_PATH', '/usr/local/bin/'); % for MAC

% make the setup files and update the input structure with the sink info
in = make_model_setup(in);  

% delete existing output files
if exist('DISC_TEST.TXT')
    delete('DISC_TEST.TXT')
    delete('DISC_DIAMETERS.TXT')
end



try        
% this runs the model 
% input files: MODEL_SETUP.TXT, SINKDIST.TXT,
fprintf('Running binary %s ...',in.binary_name);
eval(['!./' in.binary_name]);
fprintf('.done.\n');

% Saving
sav = sprintf('%s.mat',in.runfilename);
fprintf('Saving output data in %s ...', sav)
        a = load('DISC_TEST.TXT');
        b = load('DISC_DIAMETERS.TXT');
        
        out.time  = a(:,1);
        out.concs = a(:,2:end); % note: all under nucsize is rubbish
        out.drydiam = b(1,:);
        out.wetdiam = b(2,:);
        
        sav = sprintf('%s.mat',in.runfilename);
        save(sav,'in','out')
        fprintf('.done.\n');
        
catch
    fprintf('Run %s failed! error = %s',in.runfilename, lasterr)
end
    
end

%  plot_disc_auto_exact(in,out)