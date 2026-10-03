function method = selectedMethod(new)
%SELECTEDMETHOD Session state behind ptt.estimators.use and ptt.estimators.run.
%
% Private: only the functions in +ptt/+estimators can read or set it.
persistent sel
if isempty(sel), sel = ''; end
if nargin, sel = char(new); end
method = sel;
end
