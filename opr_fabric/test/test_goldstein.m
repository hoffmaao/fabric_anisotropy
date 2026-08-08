%TEST_GOLDSTEIN ptt.goldsteinFilter does what the unwrapping step needs.
%
% The NEGIS season is unwrapped by SNAPHU, and it only unwraps in reasonable
% time because run_negis_fabric.m Goldstein-filters the complex
% interferogram first. What the filter has to buy is a drop in RESIDUES -
% the 2x2 loops whose wrapped phase differences do not sum to zero - because
% those are what SNAPHU pays for in branch cuts. A filter that ran without
% reducing them would leave the unwrap as slow as it was before tiling and
% filtering, and nothing downstream would say so: the run would simply take
% an hour and be blamed on the grid size.
%
% So the test measures residues rather than asserting on the waveform, and
% pins the four properties the caller depends on:
%   1. residues fall sharply on a noisy fringe pattern, and the recovered
%      phase moves TOWARD the truth rather than away from it (a filter can
%      always cut residues by flattening the signal; that is not a pass)
%   2. alpha = 0 is exactly the identity, so the overlap-add and the taper
%      normalization reconstruct the input rather than approximating it
%   3. size and class survive (single in, single out) - the caller hands the
%      result straight to a float32 file for SNAPHU
%   4. a grid smaller than one window warns and returns the input, instead
%      of returning zeros that would read as a coherence collapse several
%      stages downstream
%
% Run: matlab -batch "run('opr_fabric/test/test_goldstein.m')"
clear;
rng(31);
t0 = tic;

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..', '..'));

Nt = 320; Nx = 256;
[cc, rr] = meshgrid(1:Nx, 1:Nt);
% A fringe pattern of the same order as the season's: a few tens of fringes
% across the grid, steeper in range than in azimuth.
phi = 2*pi * (rr/28 + cc/95);
amp = 0.6 + 0.4*exp(-rr/400);
sig = amp .* exp(1i*phi);
noise = 0.8 * (randn(Nt,Nx) + 1i*randn(Nt,Nx)) / sqrt(2);
C = single(sig + noise);

fails = 0;

% 1. residues, and phase fidelity
r_in = H_residues(angle(C));
F = ptt.goldsteinFilter(C, 0.8);
r_out = H_residues(angle(double(F)));
e_in = median(abs(angle(exp(1i*(angle(C) - phi)))), 'all');
e_out = median(abs(angle(exp(1i*(angle(double(F)) - phi)))), 'all');
drop = 1 - r_out/max(r_in, 1);
% r_in guards the test itself: a synthetic clean enough to have few
% residues to begin with would pass the drop test for free.
ok1 = r_in > 200 && drop > 0.5 && e_out < e_in;
fails = fails + ~ok1;
fprintf('residues %d -> %d (%.0f%% fewer) | phase err %.3f -> %.3f rad  %s\n', ...
  r_in, r_out, 100*drop, e_in, e_out, H_tick(ok1));

% 2. alpha = 0 is the identity
Z = ptt.goldsteinFilter(C, 0);
e_id = max(abs(double(Z(:)) - double(C(:)))) / max(abs(double(C(:))));
ok2 = e_id < 1e-9;
fails = fails + ~ok2;
fprintf('alpha 0 identity: max rel diff %.2e  %s\n', e_id, H_tick(ok2));

% 3. size and class
ok3 = isequal(size(F), size(C)) && isa(F, 'single') && ...
  isa(ptt.goldsteinFilter(double(C), 0.8), 'double') && all(isfinite(F(:)));
fails = fails + ~ok3;
fprintf('size %dx%d, class %s, all finite  %s\n', size(F,1), size(F,2), ...
  class(F), H_tick(ok3));

% 4. a grid smaller than one window is returned unfiltered, with a warning
S = C(1:40, 1:40);
ws = warning('off', 'ptt:goldsteinFilter:tooSmall');
lastwarn('');
T = ptt.goldsteinFilter(S, 0.8);
[~, wid] = lastwarn();
warning(ws);
ok4 = isequal(T, S) && strcmp(wid, 'ptt:goldsteinFilter:tooSmall');
fails = fails + ~ok4;
fprintf('40x40 grid returned unfiltered, warning "%s"  %s\n', wid, H_tick(ok4));

% 5. alpha outside [0,1] is rejected rather than silently clamped
ok5 = false;
try
  ptt.goldsteinFilter(C, 1.4);
catch err
  ok5 = strcmp(err.identifier, 'ptt:goldsteinFilter:alpha');
end
fails = fails + ~ok5;
fprintf('alpha 1.4 rejected  %s\n', H_tick(ok5));

fprintf('\n%s (%.1f s)\n', H_tick(fails == 0), toc(t0));
if fails > 0
  error('test_goldstein:failed', '%d check(s) failed', fails);
end

function n = H_residues(ph)
% Residues: the wrapped phase differences around each 2x2 loop, summed. A
% loop that does not close is a point the unwrapper has to route around.
d1 = H_wrap(ph(1:end-1, 2:end) - ph(1:end-1, 1:end-1));
d2 = H_wrap(ph(2:end,   2:end) - ph(1:end-1, 2:end));
d3 = H_wrap(ph(2:end, 1:end-1) - ph(2:end,   2:end));
d4 = H_wrap(ph(1:end-1, 1:end-1) - ph(2:end, 1:end-1));
n = nnz(abs(round((d1 + d2 + d3 + d4) / (2*pi))) > 0);
end

function w = H_wrap(x)
w = angle(exp(1i*x));
end

function s = H_tick(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end
