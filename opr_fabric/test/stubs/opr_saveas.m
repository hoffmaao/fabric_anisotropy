function opr_saveas(h_fig, fn)
% Stub of OPR opr_saveas for standalone testing
[fn_dir,~,ext] = fileparts(fn);
if ~exist(fn_dir,'dir'), mkdir(fn_dir); end
if strcmpi(ext,'.jpg') || strcmpi(ext,'.png')
  % print's JPEG device is -djpeg; -djpg is an illegal device, so every
  % .jpg the tests asked for used to fail into the warning below
  device = regexprep(lower(ext(2:end)),'^jpg$','jpeg');
  try
    print(h_fig, fn, ['-d' device], '-r100');
  catch ME
    warning('opr_saveas stub: print failed (%s)', ME.message);
  end
else
  saveas(h_fig, fn);
end
end
