# The node method

What each function in `+ptt` does and why the method is built the way it
is. Each function carries a full header of its own (`help ptt.<name>`);
this page is the map and the reasoning that connects them.
[`cresis_tutorial.md`](cresis_tutorial.md) runs it end to end, and
[`scripts.md`](scripts.md) lists the scripts around it.

## The physics in one paragraph

Ice with a horizontal crystal fabric is birefringent. Along the fabric's
slow axis a wave travels slower than across it, so the two echoes of one
reflector arrive apart by a two-way delay

    tau(z) = k * integral of dlam dz,   k = deps / (sqrt(eps_bar) c) = 0.0637 ns per metre

with dlam the difference between the two horizontal eigenvalues of the
fabric (Rathmann 2026, at nadir; the constants are in `ptt.constants`).
At Ridge A the delay reaches about 5 ns at 1500 m, almost four fringes of
the 750 MHz carrier. Everything below is a way of reading tau(z), and the
axis it is defined along, from four channels of speckle.

## Conventions

- **Axes.** `psi0` is the slow axis (the larger eigenvalue) in the antenna
  frame, clockwise, in `ptt.rotatePolarization`'s sense. The geographic
  axis is the block's heading plus `psi0`, clockwise from north, mod 180.
- **dlam.** `dlam(z) = lam(psi0) - lam(psi0 + 90)`: positive where the
  named axis is the slow one. It is not constrained to be positive.
- **Units.** Depth in metres, delay in nanoseconds, angles in radians
  inside `+ptt` and degrees in products and logs.
- **Sign of the stored data.** OPR stores the engineering sign: the phase
  of `<a conj(b)>` decreases as `a` is delayed. `ptt.xcorrTraveltime`
  conjugates into the physics sign, where it increases.
- **Angles are never unwrapped or averaged as numbers.** An axis is
  smoothed and averaged on its doubled-angle phasor.

## The functions

- `ptt.coregisterChannels` - one range delay per channel, fitted over the
  firn rows and applied to the whole record as an exact band-limited
  shift. A quad-pol set is solved through its chains (vv = V-receive +
  V-transmit, hv = V-receive, vh = V-transmit), from the HH-VV pair in the
  firn and the HV-VH pair at every depth. It is one delay and not a field
  of offsets on purpose: a coregistration that varies with depth also
  follows the delay between the slow and fast echoes and removes it. On
  Ridge A a tile search puts the median VV offset at +0.67 samples on
  north-south frames and -0.24 on east-west ones, where an instrument
  delay is one number.
- `ptt.channelGains` - the complex gains of the V receive and V transmit
  chains against H, from two constraints on the block itself: reciprocity
  (HV and VH are the same echo, at every depth) and the firn (nearly
  isotropic, so VV = HH there). Uncalibrated, the pair of polarisations
  synthesised at an azimuth is impure and the strength picks up a bias
  that alternates with the node ladder (20-40 % on a synthetic with this
  radar's gains).
- `ptt.rotatePolarization` - the channels at a synthetic antenna azimuth.
- `ptt.xcorrTraveltime` - the delay between two polarisations by windowed
  cross-correlation (50 m windows every 25 m): the envelope peak locates
  it, the carrier phase at the peak refines it, and each window also
  reports its coherence and the band's centre frequency. This is the
  delay measurement of the NEGIS two-step scheme, after Zeising et al.
  (2023).
- `ptt.crosspolNodes` - the candidate half-wave nodes. 45 deg off the
  fabric axes the cross-polarised power goes as sin^2(delta/2), so each
  extremum of the cross/co ratio against depth sits at a known delay, n
  half waves, with no tracking. About 150 m apart at Ridge A, swinging
  2-6 dB.
- `ptt.constants` - the permittivities and the speed of light behind `k`.
- `ptt.nodeFabric` - the estimator, on one along-track block:
  - **Orientation.** The strength read from the phase of the synthesised
    pair is not a cosine of the synthesis azimuth but a plateau: full
    strength within 45 deg of the slow axis, its negative within 45 deg of
    the fast axis. Fitting that square wave over 0-180 deg gives the slow
    axis directly; `plateau_ratio` is its quality.
  - **Strength.** One linear Gaussian solve for dlam(z) from the node
    delays (absolute, the fringe count) and the wrapped phase increments
    between neighbouring windows (the slope), both referenced to the firn.
    A node is read from power at zero lag, so it counts the delay in half
    waves of the band's own centre, measured per window and about 13 MHz
    below the 750 MHz LO the carrier phase counts in; the 2 % matters.
  - **The depth model.** A window's phase is the birefringent delay plus a
    phase the reflectors put between the two polarisations. That part does
    not accumulate (about 3 deg of carrier at Ridge A), so a solve that
    treats every increment as a delay hands it back as 10 m stripes. The
    model is a smooth dlam(z) (Matern 3/2, sd 0.04, 500 m, about a free
    level) plus a bounded window phase whose size each block's evidence
    sets. `resolution_m` and `n_eff` say what that leaves resolved (about
    130 m, a dozen numbers per 1.4 km column), and `intervals` gives the
    answer with no model at all.
  - **Abstention.** Below the first node where the count and the phase
    disagree, dlam is NaN, and a block with no datum at all returns NaN
    and not the prior.

## What it leaves out

- **An axis that turns with depth.** One axis per block.
- **Reflection anisotropy.** Power is used only to find the node depths.
  At Ridge A the reflection strength differs between the two axes by 2 to
  6 dB; the method does not estimate it.
- **Instrument cross-talk.** Two gains are fitted, no leakage. After them
  the cross-pol ratio still swings 0.5 dB with depth in the fabric-axis
  pair, where a clean system shows none (1.3 dB at the 45-deg pair). It
  shows as an offset of about 0.01 ns between odd and even nodes, +-0.001
  in dlam.
- **Depolarised backscatter.** Below the firn the cross-polarised power
  in the fabric's own frame is about -4 dB of the co-polarised, and nearly
  all of it is incoherent with the co-polarised echo. It lowers the
  coherence the delay is read with, which the errors carry; it is not
  modelled.
- **Ice much less coherent than Ridge A.** A window whose pair coherence
  is under 0.2 is left out (`gamma_min`). At Ridge A the coherence below
  the firn is 0.40 (quartiles 0.34 and 0.51) and 1 % of windows fall
  under the cut. On a synthetic column at 0.22, where a third of the
  windows do, the axis came back up to 15 deg off with nothing in the
  output to say so. Look at `windows.gamma` before trusting a new site.
