function opr_save(fn, varargin)
% Stub of OPR opr_save for standalone testing: gathers named variables from
% the caller workspace and saves them
s = struct();
for k = 1:numel(varargin)
  s.(varargin{k}) = evalin('caller', varargin{k});
end
[fn_dir,~,~] = fileparts(fn);
if ~exist(fn_dir,'dir'), mkdir(fn_dir); end
% -v7 so scipy/Octave/MATLAB can all read it (Octave's default is not MAT)
save(fn, '-v7', '-struct', 's');
end
