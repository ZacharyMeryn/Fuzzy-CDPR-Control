# Fuzzy Logic Control of Cable-Driven Parallel Robots in Variable Gravity Environments

Simulation code and results for the accompanying chapter, link to come upon publication



## Requirements
MATLAB -> R2018b or later 
| **Optimization Toolbox** | get 'quadprog', everything else is locally made. self check should catch you if its not properly installed



## Run order


```matlab
cd code
run_all            % runs everything in the correct order
```

or manually:

```matlab
cdpr_lowgrav_v4_4_2        % 1. main study, 15 min or so 
cdpr_dynamic_disturbance   % 2. unmodelled time-varying disturbances
cdpr_fou_noise_study       % 3. footprint width under measurement noise
```

All outputs run to 'pwd'

---

## Files

### `code/`

| file | what it does |
|---|---|
| `cdpr_lowgrav_v4_4_2.m` | Main study, contains plant, controllers, tuning, certification battery, evaluation across four gravity environments, FOU sweep, control-effort analysis, and all figures/tables |
| `cdpr_dynamic_disturbance.m` | Evaluates controllers against new time-varying disturbances that use Stribeck friction, extrusion reaction, and gust loading for the Mars setting. The goal was to create accurate depictions based on the chosen environment. no return here |
| `cdpr_fou_noise_study.m` | Sweeps the footprint half-width under measurement noise, with paired Monte Carlo duplicates |
| `run_all.m` | handles run order and prereqs |

### `results/`

`data/` holds the CSV and LaTeX exports; `figures/` holds 600 DPI PNGs. 

---

## Method summary

Four controllers are compared: an interval type-2 fuzzy controller (IT2-FLS), a
matched type-1 ablation obtained by collapsing the footprint of uncertainty to
zero width, a pole-placement PID, and a boundary-layer sliding mode controller.

Every controller outputs a desired planar force and nothing else. Downstream
they share one plant, one quadratic-program tension allocation, one saturation
limit, one velocity governor, and one gain-selection procedure, that procedure being 
a dual-plant cost, an oscillation-margin cap certified at six workspace poses, and 
a matched integral term. The point of the design is that observed differences are
attributable to the control laws rather than to unequal tuning effort.
## License

MIT — see `LICENSE`.
