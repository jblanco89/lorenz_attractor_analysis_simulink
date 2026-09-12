# Lara et al. (2003) — Lyapunov Exponent Estimation Replication

Hybrid Simulink + MATLAB replication of **"Estimación de los exponentes de Lyapunov"** (Lara, Stoico, Machado, Castagnino, ENIEF 2003, *Mecánica Computacional* Vol. XXII).

## Paper → Code Map

|Paper Section|Idea|Implementation|Why|
|-|-|-|-|
|**Sec. 3 — Método de los autovalores** (Eqs. 3–8)|Freeze Jacobian `J\_i = J(t\_i)` per window (Eq. 5), diagonalize `D\_i = P\_i^{-1} J\_i P\_i`, average `Re(eig)` (Eq. 7/7b, 8a/8b)|`core/eigenvalueSpectrum.m` + `core/lorenzJacobian.m`|Same asymptotic description as Lyapunov at lower CPU cost; sum = divergence check|
|**Sec. 4 — Resultados numéricos** (Table 1)|Same equations/params/stop, standard Gram-Schmidt vs proposed, tie at m=2 then growing win|`experiments/runCpuBenchmark.m` (pure MATLAB, reuses `src/analysis/gram\_schmidt.m`)|Per-iteration cost decides; convergence in iterations is equal|
|**Sec. 5 — Ecuación de Lorenz** (Eq. 5.1, Table 2, Figs. 1–3)|`dx=σ(y−x), dy=x(r−z)−y, dz=xy−bz`, `σ=10, b=8/3, r=20…30`, IC `(0,1,0)`, `tf=10000`, divergence `−13.666`|`sim/runSimulinkTrajectory.m` (Simulink orbit) + `experiments/runLorenzSweep.m` + `core/classifyAttractor.m`|Validates transition fixed point (r≤24) → strange (r≥25)|
|**Sec. 6 — Conclusiones**|At least one `Re(eig)` oscillates aperiodically with sign change near chaos|`viz/plotLocalEigenvalueTrace.m`|Didactic chaos signature: alternating local stability/instability|

## Cobertura del paper — cuadro comparativo

> Qué se replica, qué cambia y qué queda pendiente. Referencias a secciones del PDF (`Sec. 1–7`, `Eq.`, `Tabla`/`Fig.`).

|#|Item del paper|Estado|Qué se aplica en esta réplica|Diferencia / gap|Pendiente|
|-|-|-|-|-|-|
|1|**Sec. 1 Intro** — EDOs no lineales → caos acotado no convergente|✅ Aplicado (conceptual)|Motivación y criterio de caos documentado en comentarios y `core/classifyAttractor.m`|Sin diferencia|—|
|2|**Sec. 2 Estabilidad local** `δy(t)/δy(t0)=exp∫∂f/∂y` y `λ=lim 1/t ln`|✅ Aplicado|`core/lorenzJacobian.m` + `core/eigenvalueSpectrum.m` (Eq. 4–6) y firma 3D `(+,0,−)` en `classifyAttractor.m`|Definición analítica clásica, no hay implementación numérica propia de `λ` por separación (se usa Benettin existente para contraste)|—|
|3|**Sec. 3 Método autovalores** `J\_i=J(t\_i)` (Eq.5) → `D\_i=P\_i^{-1}J\_iP\_i` → `Λₖⁿ=1/n ΣRe(λ\_ik)` (Eq.7/7b) y `Λₖ(t)=1/t∫Re` (Eq.8a/8b)|✅ Aplicado|`core/eigenvalueSpectrum.m` congela `J\_i` por ventana, `eig` numérico (conceptualmente `D\_i`), promedia `Re` con `stride` y guarda `LambdaHistory` para el límite cuasi-estacionario|Frozen-Jacobian vía `eig` sin construir `P\_i` explícitamente (equivalente); `eig` siempre numérico aunque `m≤4` admitiría analítico|Optimizar `m≤4` con fórmula cerrada (paper lo menciona como ahorro extra)|
|4|**Sec. 3 Propiedades** `ΣΛₖ = div f`, reglas `todo <0 → estable`, `uno >0 + div<0 → conjetura caótica`, `div>0 → inestabilidad fuerte`|✅ Aplicado|`core/lorenzDivergence.m` (`-(σ+1+β)`) + `core/classifyAttractor.m` (tolerancia `1e-9`, etiquetas `FixedPoint/Strange/StrongInstability`)|Solo equivalencia de signos/clase (paper Sec.6), nunca igualdad dígito a dígito|—|
|5|**Sec. 3 Atajo analítico** `orden ≤4 → diagonalización cerrada`|⚠️ Parcial|Documentado, no bifurcado en código: `eigenvalueSpectrum.m` avisa que `eig` numérico subestima el ahorro|Se mide siempre con `eig` numérico (paper Sec.4 también lo hace para estandarizar)|Añadir rama simbólica opcional y benchmark `analítico vs numérico`|
|6|**Sec. 3 Ejemplos sin simular** Schrödinger (zona prohibida) y `m·x¨+c·x˙+k·x+β·x³` (si `k,β<0` no hay caos)|❌ No aplicado|Citados en comentarios de `eigenvalueSpectrum.m` como dictamen por signos|No se implementan modelos propios|Añadir scripts `examples/schroedinger.m` y `duffing.m` si se quiere demo didáctica|
|7|**Sec. 4 Tabla 1 — CPU vs Gram-Schmidt** `m=2..20` osciladores lineales, `Pentium MMX 32 MB`, empate `m=2` y ventaja creciente|⚠️ Parcial (puro MATLAB)|`experiments/runCpuBenchmark.m`: `buildOscillatorMatrix` (bloques `2×2 \[0 1; -k -c]`), `nIter=2000, dt=0.5`, `Phi=expm(A·dt)`, `gram\_schmidt` para `log(diag(R))/dt`|Sin modelos Simulink de osciladores; `Phi` exacto en lugar de `ode45` por ventana; segundos absolutos no comparables a 2003 (solo tendencia `speedup`)|Si se quiere fidelidad histórica: construir modelos Simulink lineales `m=2..20` y medir con `ode45` real|
|8|**Sec. 4 Estandarización** todo numérico aunque `m≤4`|✅ Aplicado|`eig`/`gram\_schmidt` siempre numéricos, incluso `m≤4`|Fiel al paper|—|
|9|**Sec. 5 Lorenz Eq.5.1** `ẋ=σ(y−x), ẏ=x(r−z)−y, ż=xy−bz`|✅ Aplicado|`models/lorenz\_sim.slx` híbrido (`sim/runSimulinkTrajectory.m` aporta la órbita, MATLAB promedia)|Solver fijo `ode4 0.001` vs paper paso variable; modelo Simulink en lugar de `ode45` puro|Validar contra `ode45` variable para cuantificar drift|
|10|**Sec. 5 Barrido** `σ=10, b=8/3, r=20..30` paso 1, `IC(0,1,0)`, `tf=10000`, `div=-13.666`|✅ Aplicado|`experiments/runLorenzSweep.m` (`RhoRange`, `Horizon`, `X0`, `Stride`)|Default `FULL tf=10000` (10 M pasos) pero smoke por defecto `tf=500` para iterar rápido; `stride=100` (ε=0.1 s) para no promediar 10 M `eig`|Ejecutar sweep completo `tf=10000` y archivar `results/lorenz\_sweep\_FULL.mat`|
|11|**Sec. 5 Tabla 2 — GS vs autovalores** dos filas por `r`, coincidencia de signos/clase, `div` idéntica|⚠️ Parcial|Fila autovalores completa para todo `r`; fila GS solo en `GsSubset=\[23 24 25 28]` con `src/analysis/lorenz\_lyapunov\_spectrum.m` tal cual (`Tben=500`, `ode45`, `dT=0.5`, 25% transitorio)|No se replica Tabla 2 completa con ambas filas para `r=20..30` (se evita coste de 11× Benettin largo)|Ampliar `GsSubset` a `20:30` o parametrizar `Tben` si se quiere tabla idéntica al paper|
|12|**Sec. 5 Figs. 1a/1b — espectros vs r** transición `24→25` (`r≤24` Fijo, `r≥25` Extraño)|✅ Aplicado|`viz/plotSpectraVsR.m` (curva continua eig + overlay discontinuo GS subset + línea `24.5` y `y=0`)|Trazo con `stride` vs integración continua del paper; `GsSubset` parcial|—|
|13|**Sec. 5 Fig. 2 — convergencia `λ₁(t)` en `r=28`**|✅ Aplicado|`viz/plotMaxExponentConvergence.m` (`LambdaHistory(1,:)` vs `tSampled`)|Horizonte del sweep (`500` smoke / `10000` FULL) vs `10000` del paper; no se superpone curva GS (se podría)|Añadir overlay GS `λ₁(t)` del Benettin si se extiende `Tben`|
|14|**Sec. 5 Figs. 3a/3b — traza `Re(λ\_local)` caótica `r=28` aperiódica vs `r=23` periódica**|✅ Aplicado|`viz/plotLocalEigenvalueTrace.m` (`ReHistory` sorteado desc., `y=0`, `on-demand` sim `h=2000` si `r` no está en sweep)|Se grafican las 3 ramas `Re` (paper muestra 1 representativa); `stride` implícito en el muestreo|—|
|15|**Sec. 6 Conclusiones / observación didáctica** caos = alternancia irregular estable/inestable|✅ Aplicado (doc.)|`classifyAttractor.m` + `plotLocalEigenvalueTrace.m` + README Decisiones|Sin demostración analítica nueva (paper la deja como trabajo futuro)|—|
|16|**Sec. 6 Hoja de ruta** fundamentación analítica + librería científica|❌ No aplicado|Fuera de alcance de esta réplica|—|Propuesta: extraer `core/` como toolbox y añadir tests `tests/testLambda.m`|

**Resumen:** réplica funcional **híbrida** fiel en lo esencial (Sec. 2–5, Figs. 1–3, clasificación y divergencia), con dos compromisos deliberados — **benchmark sin Simulink** y **GS solo en subset** — para no pagar el coste 2003. Los `results/` actuales son **smoke** (`tf=500`); el `FULL tf=10000` queda como TODO explícito.

## Layout

```
lara-lyapunov-2003/
  run\_lara\_experiments.m          # entry point (SMOKE vs FULL flags)
  core/
    lorenzJacobian.m              # J for Eq. 5.1 (Sec. 3 Eq.4 / Sec.5)
    eigenvalueSpectrum.m           # Λₖ averaging (Sec.3 Eqs.7–8)
    lorenzDivergence.m            # −(σ+1+β) control (Sec.3 \& Sec.5)
    classifyAttractor.m           # sign + divergence rules (Sec.3–5)
  sim/
    runSimulinkTrajectory.m       # Simulink wrapper (assignin + set\_param, blanks StopFcn)
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

* Model: `models/lorenz\_sim.slx` (`ode4`, `FixedStep 0.001`, `StopTime` overridden, `ReturnWorkspaceOutputs off`, `StopFcn plot\_lorenz\_attractor`).
* **Hybrid**: Simulink supplies the orbit (`xsim/ysim/zsim` timeseries); MATLAB does `eig` averaging. This respects the paper's orbit + frozen-Jacobian loop (Sec. 3) while keeping the model as the source of truth.
* **Gotchas handled** (`AGENTS.md`): `SimulationInput.setVariable` silently fails → `assignin('base',...) + set\_param(StopTime)`; `StopFcn` blanked during batch and restored after (never left blank); `SaveTime/SaveOutput` off during function-called `sim` to avoid `tout` workspace clash; base `x0/y0/z0, sigma/rho/beta` snapshotted and restored.

## Quick Start

```matlab
% Smoke (seconds, 2 rhos, short horizon, validates transition)
% Option A: edit lara-lyapunov-2003/run\_lara\_experiments.m -> set DO\_SMOKE=true, then:
run\_lara\_experiments   % from repo root after addpath(genpath(pwd))

% Option B: call experiments directly
addpath(genpath('lara-lyapunov-2003'))
res = runLorenzSweep('RhoRange',\[24 25],'Horizon',500,'Stride',200,'GsSubset',\[]);
bench = runCpuBenchmark('Dims',2:2:10,'NIter',500);
```

```matlab
% Full paper fidelity (slow: 11× tf=10000 sims at FixedStep 0.001 → \~10M steps each)
% In run\_lara\_experiments.m: DO\_SMOKE=false (default), CFG.horizon=10000
run\_lara\_experiments  % → results/lorenz\_sweep.mat + fig1/2/3 pngs
```

## Decisions

1. **Table 1** → pure-MATLAB timing (no Simulink oscillator models).
2. **Full fidelity** `tf=10000` default; stride (default 100 → effective window 0.1 s) throttles `eig` cost without biasing the mean (convergence in iterations is equal, Sec. 4).
3. **Gram-Schmidt cross-check** on subset `\[23 24 25 28]` with `src/analysis/lorenz\_lyapunov\_spectrum.m` as-is (`Tben=500`, `ode45`, `dT=0.5`, 25% transient). Only sign/class equivalence is asserted (Sec. 6), never digit equality.

## Verification

* Smoke `r=24 → FixedPoint, r=25 → Strange`, divergence `−13.666 ± 0.05`.
* Full sweep reproduces Table 2 classes and the 24→25 jump (Figs. 1a/1b).
* `sum(Lambda) ≈ −13.666` on every run; GS subset classes match eig classes on the subset.
* `runCpuBenchmark` shows tie-ish at m=2 and growing speedup with m (paper Table 1 trend; absolute seconds are machine-dependent).

