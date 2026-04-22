# VertexModelChiralWetEdges

A 2D vertex model of epithelial tissues with **chiral activity** and **wet-edge dynamics**. Built on the vertex model framework originally written by Matej Krajnc et al. (please acknowledge in derived work).

## Repository layout

| Path | Contents |
| --- | --- |
| [main.cpp](main.cpp) | Entry point. Reads simulation parameters from the command line, initializes the tissue, runs the time loop, and writes outputs. |
| [vertX2D/](vertX2D/) | Core C++ headers: forces (area, perimeter, chiral, nematic), T1 transitions, equation of motion, initial conditions, output, geometry utilities. |
| [seeds/](seeds/) | 100 pre-generated disorder seeds (`seed_1.dat` … `seed_100.dat`) used to initialize tissues. |
| [Python_code/](Python_code/) | Post-processing and analysis: Kirkwood stress, velocity averages (`V`, `Vr`, `Vphi`), vertex density, MSD, azimuthal averaging, ensemble analysis, plotting. |
| `main.opp` | Compiled binary (built from `main.cpp`). |

## Build

Requires a C++ compiler with C++17 (`std::filesystem`) and the **Eigen** sparse linear algebra library (used for the wet-dynamics connectivity matrix).

```bash
g++ -std=c++17 -O3 -I/path/to/eigen main.cpp -o main.opp
```

## Run

The binary expects **23 command-line arguments** in the following order:

| # | Name | Description |
| --- | --- | --- |
| 1 | `Nx` | Initial system size in cells along x (12 for quick tests, 32 for production). |
| 2 | `c_kAb` | Target-area modulus (e.g. `0.5`). |
| 3 | `kPer` | Target-perimeter modulus (e.g. `0.01`). |
| 4 | `P0` | Target perimeter (typically `3.6`–`3.95`). |
| 5 | `disorder` | Disorder amplitude (`0.001` = hexagonal, `0.02` = moderate, `0.03` = high). |
| 6 | `tMAX` | Total simulation time (e.g. `1000`). |
| 7 | `delta_write_time` | Interval between scalar-output writes (e.g. `10`). |
| 8 | `delta_image_time` | Interval between full-tissue snapshots (e.g. `100` or `1000`). |
| 9 | `chiral_alpha_1` | Chiral activity, population 1 (inside the droplet). |
| 10 | `chiral_alpha_2` | Chiral activity, population 2. |
| 11 | `initial_condition` | `0`: drop, `1`: random mix, `2`: split along x, `3`: split along y, `4`: chick (streak + node). |
| 12 | `initial_parameter` | IC parameter: drop radius (fraction of half-box-y) or population-2 fraction for mix. |
| 13 | `initial_seed` | RNG seed for the initial condition. |
| 14 | `tiling_seed` | RNG seed for the disordered tiling. |
| 15 | `dynamics_type` | `0`: dry, `1`: vertex–vertex wet, `2`: edge-wet. |
| 16 | `_gamma_wet` | External friction (must be non-zero). |
| 17 | `_zeta_wet` | Internal friction (often `1`, with `_gamma_wet` one or more orders of magnitude smaller). |
| 18 | `transition_type` | Drop interface: `0`: sharp, `1`: smooth. |
| 19 | `streak_length` | Chick IC: streak length. |
| 20 | `streak_width` | Chick IC: streak width. |
| 21 | `streak_center_y` | Chick IC: streak y-center. |
| 22 | `node_radius` | Chick IC: node radius. |
| 23 | `node_center_y` | Chick IC: node y-center. |

Example:

```bash
./main.opp 32 0.5 0.01 3.85 0.02 1000 10 100 \
           0.1 0.0 0 0.5 1 1 2 0.1 1.0 0 \
           0 0 0 0 0
```

## Output

Each run creates a directory `./output/out_Nx_<Nx>_kA_<kA>_kP_<kP>_P0_<P0>_…/` encoding all parameters in the name, containing:

- `energies.dat` — energy time series
- `sorting.dat` — population-sorting metric
- `msd.dat` — mean squared displacement
- `max_move.dat` — maximum per-step vertex displacement
- `chiral_area.dat` — chiral-area time series
- `T1_transitions.dat` — T1 event log
- Periodic full-tissue snapshots (vertices, cells, velocities, wet state, etc.)

## Analysis (Python)

The [Python_code/](Python_code/) directory contains readers and analysis utilities that consume the per-run output folders. Highlights:

- **Readers**: `read_file_folder*.py` parse snapshots for vertices, velocities, cell spins, wet state, T1 events, etc.
- **Fields**: `V_ave_func.py`, `Vr_ave_func.py`, `Vphi_ave_func.py`, `LapV_ave_func.py`, `vertex_density_ave_func.py` compute time/azimuthal averages.
- **Stress**: `Kirkwood_stress*.py` and `Kirkwood_stress.ipynb` compute Kirkwood stress and dissipation, including gradient and time-average plots.
- **Ensemble**: `ensemble_average_analysis.py` aggregates across seeds.
- **Fits / plots**: `fit_vphi.py`, `plot_vphi.py`, `numerical_differentiation.py`, `wet_force_balance_out/` (generated figures).

## Notes

- Global defaults at the top of [main.cpp](main.cpp) are left unchanged by design — simulation parameters flow in through argv.
- `dynamics_type == 1` builds a sparse connectivity matrix via `update_M()` (Eigen `SparseLU`); expect higher memory/CPU than dry runs.
- The `*.sh` driver scripts and `output/` directories are intentionally gitignored.
