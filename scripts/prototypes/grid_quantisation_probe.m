%GRID_QUANTISATION_PROBE Does the variable-projection grid limit accuracy?
%
% ptt.quadpolFabricLS is a separable nonlinear least squares solved by
% variable projection: three non-linear parameters (theta_0, delta_0,
% delta') on an exhaustive grid, four linear ones eliminated in closed
% form at every node. The grid is 60 x 26 x 24 = 37440 four-by-four solves
% per window, and delta_0 is quantised at 15 degrees, which looks coarse.
%
% IT IS NOT COSTING ANYTHING. MEASURED, refining delta_0 from 15 to 5 to
% 2.5 degrees - six times the nodes - moves neither the axis nor the
% contrast at all:
%
%   as shipped   (th 3, d0 15)    theta err 0.015 deg   dlam 0.0499
%   finer d0     (th 3, d0 5)     theta err 0.015 deg   dlam 0.0499
%   finer d0     (th 3, d0 2.5)   theta err 0.015 deg   dlam 0.0499
%                                                       truth  0.0500
%
% An axis error of 0.015 degrees on a 3 degree grid also says the theta
% grid is not the limit: the cost surface is smooth enough that the grid
% minimum lands essentially on the true one.
%
% WHY THIS SCRIPT EXISTS. Optimiser improvements are the tempting kind of
% work - Golub-Pereyra's analytic Jacobian with Levenberg-Marquardt,
% Kaufman's simplification, or eliminating delta_0 analytically (the model
% IS linear in cos delta_0 and sin delta_0 once the denominator is taken
% from the data, which would collapse the grid 24-fold). All are real
% techniques and none of them would improve a number this project reports,
% because the solver is already exact where its model holds. Run this
% before spending time there.
%
% What DOES limit the answers, all measured elsewhere in this repo:
% decorrelation (at EastGRIP dlam tracks coherence at +0.54 down the
% column), the constant-eigenframe assumption (6.9x misfit on a column
% whose axis turns 60 degrees; the once-reported 42.6 degree gap between
% held and free axes was a units error and is 1.7 degrees), the antenna
% pedestal (which antenna-locks every
% power-only route on this system), and a linearised posterior that runs
% about 2x optimistic.
%
if ~exist('ptt_root','var'), ptt_root = fileparts(fileparts(fileparts(mfilename('fullpath')))); end % this checkout's root, holding +ptt
addpath(ptt_root);
rng(5);
fc=750e6; z=(10:10:1200).'; Nx=400; DL=0.05; TH=deg2rad(35);
C=ptt.constants(); n=sqrt(C.eps_bar);
gpd=2*pi*fc*2*(C.deps/(2*n))/(C.c*1e9);
lay=struct('top_m',0,'dlam',DL,'theta',TH,'r_db',0,'gx_db',0);
fm=ptt.fujitaModel(lay,z,0,struct('fc',fc,'win_m',30));
r=(randn(numel(z),Nx)+1i*randn(numel(z),Nx))/sqrt(2);
S=struct('hh',fm.s_hh(:,1).*r,'vv',fm.s_vv(:,1).*r, ...
         'hv',fm.s_hv(:,1).*r,'vh',fm.s_hv(:,1).*r);
base=struct('fc',fc,'psi_step_deg',3,'win_fit_m',120,'step_m',60, ...
  'deramped',true,'pedestal','window');
% expected axis is TH + 90 (the convention pinned in test_ls_vs_fujita)
TH_EXP = TH + pi/2;
fprintf('%-34s %9s %9s %9s\n','grid','th err','dlam','nodes');
cfg = { 'as shipped (th 3, d0 15)',      3, 24, 26; ...
        'finer d0 only (th 3, d0 5)',    3, 72, 26; ...
        'finer d0 (th 3, d0 2.5)',       3, 144, 26 };
for k=1:size(cfg,1)
  o=base; o.theta_step_deg=cfg{k,2};
  o.d0_grid=(0:(360/cfg{k,3}):360-(360/cfg{k,3}))*pi/180;
  o.dd_grid=linspace(0, 0.25*gpd, cfg{k,4});
  t0=tic; out=ptt.quadpolFabricLS(S,z,o); el=toc(t0);
  g=isfinite(out.theta0);
  th=angle(median(exp(2i*out.theta0(g))))/2;
  e=rad2deg(abs(angle(exp(2i*(th-TH_EXP)))/2));
  dl=median(out.dlam(isfinite(out.dlam)));
  fprintf('%-34s %9.3f %9.4f %9d   (%.0f s)\n', cfg{k,1}, e, dl, ...
    round(180/cfg{k,2})*cfg{k,3}*cfg{k,4}, el);
end
fprintf('\ntruth: theta %.1f deg (expected %.1f after the +90 convention), dlam %.4f\n', ...
  rad2deg(TH), mod(rad2deg(TH_EXP),180), DL);
