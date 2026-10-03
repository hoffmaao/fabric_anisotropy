function method = current()
%CURRENT The estimation method selected for this session ('' if none).
%
% method = ptt.estimators.current()
%
% See also ptt.estimators.use, ptt.estimators.run.
method = selectedMethod();
end
