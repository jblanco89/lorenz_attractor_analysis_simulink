# Lara et al. (2003) — Lyapunov Exponent Estimation Replication

Hybrid Simulink + MATLAB replication of **"Estimación de los exponentes de Lyapunov"** (Lara, Stoico, Machado, Castagnino, ENIEF 2003, *Mecánica Computacional* Vol. XXII).

## Paper → Code Map

| Paper Section | Idea | Implementation | Why |
|---|---|---|---|
| **Sec. 3 — Método de los autovalores** (Eqs. 3–8) | Freeze Jacobian `J_i = J(t_i)` per window (Eq. 5), diagonalize `D_i = P_i^{-1} J_i P_i`, average `Re(eig)` (Eq. 7/7b, 8a/8b) | `core/eigenvalueSpectrum.m` + `core/lorenzJacobian.m` | Same asymptotic description as Lyapunov at lower CPU cost; sum = divergence check |
| **Sec. 4 — Resultados numéricos** (Table 1) | Same equations/params/stop, standard Gram-Schmidt vs proposed, tie at m=2 then growing win | `experiments/runCpuBenchmark.m` (pure MATLAB, reuses `src/analysis/gram_schmidt.m`) | Per-iteration cost decides; convergence in iterations is equal |
| **Sec. 5 — Ecuación de Lorenz** (Eq. 5.1, Table 2, Figs. 1–3) | `dx=σ(y−x), dy=x(r−z)−y, dz=xy−bz`, `σ=10, b=8/3, r=20…30`, IC `(0,1,0)`, `tf=10000`, divergence `−13.666` | `sim/runSimulinkTrajectory.m` (Simulink orbit) + `experiments/runLorenzSweep.m` + `core/classifyAttractor.m` | Validates transition fixed point (r≤24) → strange (r≥25) |
| **Sec. 6 — Conclusiones** | At least one `Re(eig)` oscillates aperiodically with sign change near chaos | `viz/plotLocalEigenvalueTrace.m` | Didactic chaos signature: alternating local stability/instability |

## Layout

```
lara-lyapunov-2003/
  run_lara_experiments.m          # entry point (SMOKE vs FULL flags)
  core/
    lorenzJacobian.m              # J for Eq. 5.1 (Sec. 3 Eq.4 / Sec.5)
    eigenvalueSpectrum.m           # Λₖ averaging (Sec.3 Eqs.7–8)
    lorenzDivergence.m            # −(σ+1+β) control (Sec.3 & Sec.5)
    classifyAttractor.m           # sign + divergence rules (Sec.3–5)
  sim/
    runSimulinkTrajectory.m       # Simulink wrapper (assignin + set_param, blanks StopFcn)
  experiments/
    runLorenzSweep.m              # r=20…30 sweep → Table 2
    runCpuBenchmark.m             # m=2…20 timing → Table 1
  viz/
    plotSpectraVsR.m              # Figs. 1a/1b
    plotMaxExponentConvergence.m  # Fig. 2
    plotLocalEigenvalueTrace.m    # Figs. 3a/3b
  results/                        # .mat + exported pngs (git-ignored except README)
```

## Simulink Use

- Model: `models/lorenz_sim.slx` (`ode4`, `FixedStep 0.001`, `StopTime` overridden, `ReturnWorkspaceOutputs off`, `StopFcn plot_lorenz_attractor`).
- **Hybrid**: Simulink supplies the orbit (`xsim/ysim/zsim` timeseries); MATLAB does `eig` averaging. This respects the paper's orbit + frozen-Jacobian loop (Sec. 3) while keeping the model as the source of truth.
- **Gotchas handled** (`AGENTS.md`): `SimulationInput.setVariable` silently fails → `assignin('base',...) + set_param(StopTime)`; `StopFcn` blanked during batch and restored after (never left blank); `SaveTime/SaveOutput` off during function-called `sim` to avoid `tout` workspace clash; base `x0/y0/z0, sigma/rho/beta` snapshotted and restored.

## Quick Start

```matlab
% Smoke (seconds, 2 rhos, short horizon, validates transition)
% Option A: edit lara-lyapunov-2003/run_lara_experiments.m -> set DO_SMOKE=true, then:
run_lara_experiments   % from repo root after addpath(genpath(pwd))

% Option B: call experiments directly
addpath(genpath('lara-lyapunov-2003'))
res = runLorenzSweep('RhoRange',[24 25],'Horizon',500,'Stride',200,'GsSubset',[]);
bench = runCpuBenchmark('Dims',2:2:10,'NIter',500);
```

```matlab
% Full paper fidelity (slow: 11× tf=10000 sims at FixedStep 0.001 → ~10M steps each)
% In run_lara_experiments.m: DO_SMOKE=false (default), CFG.horizon=10000
run_lara_experiments  % → results/lorenz_sweep.mat + fig1/2/3 pngs
```

## Decisions

1. **Table 1** → pure-MATLAB timing (no Simulink oscillator models).
2. **Full fidelity** `tf=10000` default; stride (default 100 → effective window 0.1 s) throttles `eig` cost without biasing the mean (convergence in iterations is equal, Sec. 4).
3. **Gram-Schmidt cross-check** on subset `[23 24 25 28]` with `src/analysis/lorenz_lyapunov_spectrum.m` as-is (`Tben=500`, `ode45`, `dT=0.5`, 25% transient). Only sign/class equivalence is asserted (Sec. 6), never digit equality.

## Verification

- Smoke `r=24 → FixedPoint, r=25 → Strange`, divergence `−13.666 ± 0.05`.
- Full sweep reproduces Table 2 classes and the 24→25 jump (Figs. 1a/1b).
- `sum(Lambda) ≈ −13.666` on every run; GS subset classes match eig classes on the subset.
- `runCpuBenchmark` shows tie-ish at m=2 and growing speedup with m (paper Table 1 trend; absolute seconds are machine-dependent).

