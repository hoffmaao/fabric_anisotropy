function ang_deg = axisFold(theta, ref)
%AXISFOLD Angle between a fabric AXIS and a reference, folded to [0, 90].
%
% Shared by ptt.orientationToFlow and ptt.orientationToStrain so the two
% cannot drift apart on the one piece of arithmetic that is easy to get
% wrong.
%
% A fabric axis is defined modulo 180 degrees: it has no head or tail, so
% "aligned" and "anti-aligned" are the same statement. Doubling the angle
% folds those two onto the same point, and halving the wrapped difference
% lands the answer in [0, 90] - 0 parallel, 90 perpendicular - with no
% unwrapping anywhere, which is this project's standing rule for angles.
%
% The same fold is correct whether the reference is a DIRECTION (flow,
% modulo 360) or itself an AXIS (a strain-rate principal axis, modulo
% 180), because doubling collapses a direction onto its own axis. What
% differs is only the MEANING of the answer, which is the caller's to
% document, not this function's to decide.
ang_deg = abs(rad2deg(angle(exp(2i * (double(theta) - double(ref)))) / 2));
end
