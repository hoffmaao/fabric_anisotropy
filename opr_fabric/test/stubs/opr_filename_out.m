function fn = opr_filename_out(param, out_path)
% Stub of OPR opr_filename_out for standalone testing: paths under
% param.stub_out_root instead of gRadar.out_path
fn = fullfile(param.stub_out_root, ['CSARP_' out_path], param.day_seg);
end
