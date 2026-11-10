# Pedagogical OLG model — 80 annual cohorts

Ages 20–99, flexible prices, retirement at 65, one asset, endogenous capital, wage and interest rate.
Put **Dynare 7.x (at least 7.1; 7.2 recommended)** on the MATLAB path, then run from this directory.

Download the seven code files listed below into one folder and run from that folder. MAT results are saved alongside the `.mod` files, with no `results/` subfolder.

## Files

Each model is a single `.mod` file read from top to bottom: calibration, model, steady state, shock,
simulation, figures. The unconstrained `steady_state_model` block displays the analytic household life cycle directly in `olg80.mod`; a small helper solves only the scalar market-clearing rate. The constrained case uses a separate household solver.

| File | Role |
|---|---|
| `olg80.mod` | No borrowing limit: analytic life cycle, scalar market clearing, calibration figures, and analytic windfall MPCs saved to `mpc_bc0.mat` (TFP run) |
| `olg80_interest_rate.m` | Scalar asset-market root; evaluates analytic household saving at each trial rate while `olg80.mod` displays the formula for first-cohort consumption and builds the entire steady state |
| `rank80.mod` | Representative-agent benchmark; draws the OLG-vs-RANK comparison (four aggregate subplots) |
| `compare_olg_rank.m` | Student launcher: choose the shock and run OLG then RANK, or redraw their saved comparison |
| `olg80_bc.mod` | Appendix: borrowing limit `a_j >= 0` only; windfall MPCs saved to `mpc_bc1.mat` |
| `olg80_bc_steadystate.m`, `lifecycle_bc.m` | Constrained steady state and prefix-budget solver; fixed-price MPCs recompute remaining life after a current windfall (`dW_j` remains available for separate GE experiments) |

## Slide figures (separate scripts in `plots/`)

Each lecture figure has its own standalone script: it first runs the required Dynare `.mod` files, then draws and exports the PDF/PNG used in the slides. Keep each script beside the downloaded model files (or in `model/plots/`): it locates the model files itself and only needs Dynare 7.x on the MATLAB path. In a single download folder, the figures are saved in that folder (in `figures/` for the borrowing-limit script). `plot_kaplan2014_mpc.m` plots hardcoded published estimates and needs no model run.

| Script | Figure produced |
|---|---|
| `plots/plot_olg80_figures.m` | `olg80_calibration_profiles.pdf`, `olg80_lifecycle.pdf`, `olg80_asset_market.pdf` |
| `plots/plot_olg80_rank_comparison.m` | `olg80_tfp_transition.pdf`, `olg80_depreciation.pdf` |
| `plots/plot_olg80_tfp_consumption_focus.m` | `olg80_tfp_consumption_focus.png` (consumption on impact by age, TFP shock) |
| `plots/plot_olg80_borrowing_limit.m` | `olg80_borrowing_mpc.pdf`, `olg80_tfp_transition_bc.pdf` (appendix) |
| `plots/plot_kaplan2014_mpc.m` | `kaplan2014_mpc.pdf` (empirical MPC evidence) |

## Running order

Start with `compare_olg_rank.m` in MATLAB: `shock_type = 1` compares the productivity shock, `shock_type = 2` the depreciation shock. `rerun_models = true` (default) runs both models; set it to `false` to compare existing MAT files without rerunning. OLG also displays its own aggregate responses (consumption, output, capital and net interest rate) in a separate 2×2 figure.

The equivalent individual commands are:

```matlab
dynare olg80 -DSHOCK_TYPE=1 noclearall     % TFP shock
dynare rank80 -DSHOCK_TYPE=1 noclearall    % RANK
dynare olg80 -DSHOCK_TYPE=2 noclearall     % depreciation shock (financial crisis)
dynare rank80 -DSHOCK_TYPE=2 noclearall 
```

Appendix (after the TFP runs above):

```matlab
dynare olg80_bc noclearall                        % constrained MPC figure (subplots)
```

`rank80.mod` reads the OLG MAT file for the selected shock (`results_tfp.mat` or `results_dep.mat`), so run `olg80.mod` first for that shock.
Its TFP run also writes the analytic no-limit MPC profile (`mpc_bc0.mat`);
`olg80_bc.mod` then writes the constrained profile (`mpc_bc1.mat`).
Both MPC profiles hold wages and the interest rate fixed at each model's own stationary equilibrium.
The old `-DBC=0/1` switch has been removed: choose the corresponding model file instead.

## Calibration

The shared preference, production and TFP parameters in `olg80.mod`, `rank80.mod` and `olg80_bc.mod` use the same calibration. To change one, edit
its value in all three files in a copy of the folder (the constrained model has additional windfall dummies). The scalar-root helper and constrained steady-state function read the executed
`M_.params` and keep no second calibration. Equal cohort masses are assumed.
The free household profile is displayed in `olg80.mod` and follows the closed-form lifetime budget and Euler growth;
the constrained profile instead enforces every remaining prefix budget.

To inspect the macro expansion with `J=3`: `dynare olg80 onlymacro savemacro`.
