%POSTERIOR_FROM_GRID_PROBE Can the existing grid be read as a posterior?
%
% Asked while considering a Hamiltonian Monte Carlo framing of the fabric
% inversion, to get honest joint uncertainty on direction and strength.
%
% THE FINDING: for the per-window problem, HMC is the wrong tool because
% the posterior is already being computed. ptt.quadpolFabricLS evaluates a
% 60 x 26 x 24 grid over the three non-linear parameters and eliminates the
% four linear ones exactly at every node. Normalising exp(-cost/2 sigma^2)
% over that grid IS the joint posterior of (theta_0, delta_0, delta') -
% exactly, with no sampler, no step-size tuning and no convergence
% diagnostics to defend. Three parameters is squarely grid territory.
%
% THE ONE THING THAT HAD TO BE CHECKED, and this script checks it. The
% estimator minimises the PROFILE cost at the plug-in nuisances, whereas
% integrating them out adds an Occam term:
%     -2 log p(d|m_n) = cost(m_n)/sigma^2 + log det M(m_n)
% If log det M varied with theta_0, the grid would be the wrong shape and
% the profile cost could not be read as a posterior. MEASURED over the
% theta_0 grid: the profile cost spans 1.13e5 while log det M spans 0.135,
% a ratio of 1.2e-6, and the minimum sits at 36.0 deg either way. The
% Occam term is negligible here, so the grid can be normalised as it
% stands.
%
% WHAT IS STILL MISSING is the noise scale. The CRB weights are relative,
% so exp(-cost/2 sigma^2) has no absolute width until sigma^2 is supplied
% or marginalised (a scale prior gives a Student-t marginal). That, not the
% sampler, is the work standing between here and a defensible interval.
%
% WHERE HMC WOULD EARN ITS PLACE: the full-profile problem in
% ptt.fabricGLS, which carries 3 x n_layer parameters where a grid is
% impossible, and whose linearised posterior test_fabric_gls MEASURES to be
% about 2x optimistic. That is a real gap and a real use for a sampler.
%
% WHAT NO SAMPLER FIXES: a misspecified likelihood. HMC would faithfully
% sample the wrong posterior if the effective look count or the correlation
% between samples is wrong, and it does nothing about the constant-
% eigenframe assumption.
%
% theta_0, or only shift it? The estimator minimises the PROFILE cost
%   cost(m_n) = C2 - 2 x'r + x'Mx
% at the plug-in x. Integrating x out instead adds an Occam term:
%   -2 log p(d|m_n) = cost/sigma^2 + log det M(m_n)
% If log det M is flat in theta_0 the two agree up to a constant and the
% grid is already a posterior; if it is not, a Bayesian treatment must
% carry it.
if ~exist('ptt_root','var'), ptt_root = fileparts(fileparts(fileparts(mfilename('fullpath')))); end % this checkout's root, holding +ptt
addpath(ptt_root);
rng(3);
fc=750e6; z=(400:2:900).'; psi=(0:3:177)*pi/180;
TH=deg2rad(35); DL=0.05;
C=ptt.constants(); nice_=sqrt(C.eps_bar);
gpd=2*pi*fc*2*(C.deps/(2*nice_))/(C.c*1e9);
dd=gpd*DL; zc=mean(z); u=(z-zc).';
% synthetic coherence field from the same closed form the estimator fits
mk=@(th,d0) buildC(psi,u,th,d0,dd);
Cw = mk(TH, 1.1);
Cw = Cw .* (0.85 + 0.10*randn(size(Cw)));           % speckle-ish
W  = abs(Cw).^2 ./ max(1-abs(Cw).^2, 0.02);
sb = sin(2*psi(:)); s2b = sb.^2;
wrow = sum(W,2);
K.W2 = (sb.^2).'*wrow; K.W3=(sb.^3).'*wrow; K.W4=(sb.^4).'*wrow;
K.tau = 1e-3*max(K.W2,realmin);
C2 = sum(W.*abs(Cw).^2,'all');
ths = (0:2:178)*pi/180;
cost=nan(size(ths)); ld=nan(size(ths));
for i=1:numel(ths)
  H = buildC(psi,u,ths(i),1.1,dd);
  rH=real(H); iH=imag(H);
  G1=sum(W.*(real(Cw).*rH+imag(Cw).*iH),'all');
  G2=sum(W.*(rH.^2+iH.^2),'all');
  rr=sum(W.*rH,2); ii=sum(W.*iH,2);
  b01=sb.'*rr; b02=sb.'*ii; b03=s2b.'*rr;
  M=[G2 b01 b02 b03; b01 K.W2+K.tau 0 K.W3; b02 0 K.W2+K.tau 0; b03 K.W3 0 K.W4+K.tau];
  r=[G1; sb.'*sum(W.*real(Cw),2); sb.'*sum(W.*imag(Cw),2); s2b.'*sum(W.*real(Cw),2)];
  x=M\r;
  cost(i)=C2-2*(x.'*r)+x.'*M*x;
  ld(i)=log(det(M));
end
fprintf('over the theta_0 grid:\n');
fprintf('  profile cost   range %.4g  (min %.4g)\n', max(cost)-min(cost), min(cost));
fprintf('  log det M      range %.4g\n', max(ld)-min(ld));
fprintf('  ratio (logdet range / cost range) = %.3g\n', (max(ld)-min(ld))/(max(cost)-min(cost)));
[~,ia]=min(cost); [~,ib]=min(cost+ld);
fprintf('  argmin profile %.1f deg    argmin profile+logdet %.1f deg\n', ...
  rad2deg(ths(ia)), rad2deg(ths(ib)));
function H = buildC(psi,u,th,d0,dd)
mu=cos(2*(psi(:)-th)); A=(1-mu.^2)/2; B=(1+mu.^2)/2;
d=d0+dd*u; cd=cos(d); sd=sin(d);
H=(A+B.*cd+1i*(mu.*sd))./max(1-A.*(1-cd),0.05);
end
