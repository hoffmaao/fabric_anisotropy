function opr_saveas(h_fig, fn)
% Stub of OPR opr_saveas for standalone testing
[fn_dir,~,ext] = fileparts(fn);
if ~exist(fn_dir,'dir'), mkdir(fn_dir); end
if strcmpi(ext,'.jpg') || strcmpi(ext,'.png')
  try
    print(h_fig, fn, ['-d' lower(ext(2:end))], '-r100');
  catch ME
    warning('opr_saveas stub: print failed (%s)', ME.message);
  end
else
  saveas(h_fig, fn);
end
end
