# Pedagogical OLG model — 80 annual cohorts

Ages 20–99, flexible prices, retirement at 65, one asset, endogenous capital, wage and interest rate.
Put **Dynare 7.x (at least 7.1; 7.2 recommended)** on the MATLAB path, then run from this directory.

## Files

Each model is a single `.mod` file read from top to bottom: calibration, model, steady state, shock,
simulation, figures. Only the steady state needs a MATLAB function.

| File | Role |
|---|---|
| `olg80.mod` | 80-cohort model; also draws the calibration, life-cycle and asset-market figures |
| `olg80_steadystate.m` | stationary equilibrium (clears the asset market) |
| `rank80.mod` | representative-agent benchmark; draws the OLG-vs-RANK figure of the same experiment |
| `olg80_bc.mod` | appendix: borrowing limit (`-DBC=1`) and windfall experiment (`-DSHOCK_TYPE=2`) |
| `olg80_bc_steadystate.m`, `lifecycle_bc.m` | appendix steady state and the household life-cycle solver |

## Running order

```matlab
dynare olg80 -DSHOCK_TYPE=1 noclearall     % TFP shock; also draws the calibration figures
dynare rank80 -DSHOCK_TYPE=1 noclearall    % RANK, figure olg80_tfp_transition.pdf
dynare olg80 -DSHOCK_TYPE=2 noclearall     % depreciation shock (financial crisis)
dynare rank80 -DSHOCK_TYPE=2 noclearall    % figure olg80_depreciation.pdf
```

Appendix (after the TFP runs above):

```matlab
dynare olg80_bc -DBC=1 -DSHOCK_TYPE=1 noclearall                 % MPC and TFP-with-limit figures
dynare olg80_bc -DBC=0 -DSHOCK_TYPE=2 -DTARGET=6 noclearall      % windfall to age 25, no limit
dynare olg80_bc -DBC=1 -DSHOCK_TYPE=2 -DTARGET=6 noclearall      % same, with the limit
```

`TARGET` is the cohort index (age = TARGET + 19; 56 is age 75). `rank80.mod` reads the OLG steady
state from `results/results_tfp.mat`, so run `olg80.mod` first. Figures go to `../figures/`.

## Calibration

The calibration block at the top of `olg80.mod`, `rank80.mod` and `olg80_bc.mod` is identical (the
appendix file adds a transfer process in place of the depreciation process). To change a parameter, edit
the block in all three files in a copy of the folder. The steady-state functions read the executed
`M_.params` and keep no second calibration. Equal cohort masses are assumed. The
`constrained` parameter can trigger a Dynare unused-parameter warning: the steady-state function reads it,
while `BC` selects the equations at macro time.

To inspect the macro expansion with `J=3`: `dynare olg80 onlymacro savemacro`.
