# ElasticLoading

**MATLAB tools for computing elastic surface-loading deformation with PREM Green's functions**

ElasticLoading is a small research-oriented MATLAB toolkit for estimating the **vertical (Up), eastward (East), and northward (North) elastic displacements** caused by changes in surface mass loading. It uses tabulated Green's functions for the Preliminary Reference Earth Model (PREM), following a Farrell-style loading formulation.

The repository currently contains a generic pointwise load-convolution routine and a worked application that converts **global CSR GRACE mascon equivalent-water-height (EWH) trends** into regional elastic displacement rates. **The example is configured for the North Pacific (January 2005–December 2024)**, but its paths and regional polygon can be changed.

> **Research status:** The supplied scripts are intended for research and adaptation, not a turnkey or independently validated processing pipeline. External input data and some dependencies are not provided in this repository.

## Contents

| File | Purpose |
| --- | --- |
| [`cf_PREM.txt`](cf_PREM.txt) | Tabulated PREM Green's-function values. Column 1 gives angular distance; columns 4 and 5 are used for vertical and horizontal response, respectively. |
| [`mass_load_hv.m`](mass_load_hv.m) | Computes Up/East/North displacement in **mm** at one observation location from a regular-grid surface load. |
| [`green_hv_region.m`](green_hv_region.m) | Applies `mass_load_hv` to a list of observation locations. |
| [`P08_CSR_mascon_elastic_seafloor_PREM.m`](P08_CSR_mascon_elastic_seafloor_PREM.m) | Example function: fits global CSR mascon EWH trends over 2005–2024 and calculates regional elastic displacement **rates in mm/yr**, then exports MAT, CSV, GMT input, and a map. |
| [`LICENSE`](LICENSE) | MIT license. |

## Method

For a load cell `j`, the input EWH `h_j` (m) is converted to mass by

```text
ΔM_j = ρw × h_j × A_j
```

where `ρw = 1000 kg/m³` is the configured water density and `A_j` (m²) is the spherical grid-cell area. The existing implementation evaluates tabulated Green's functions at the **angular distance between the source cell centre and the observation point**, and sums the individual cell responses:

```text
u_up = Σ_j G_up(ψ_j) ΔM_j
```

The code uses the normalization

```text
G_up(ψ) = ur2(ψ) / [10^12 × R × ψ]
G_h (ψ) =  v2(ψ) / [10^12 × R × ψ]
```

with `R = 6,371,000 m` and **`ψ` in radians in the denominator**; the table lookup uses angular degrees in the default settings. East and north components are obtained using the source-to-receiver azimuth and the horizontal Green's function, with the sign conventions in the MATLAB source. The displacement outputs are ordered **[Up, East, North]**.

The CSR example independently fits a linear trend at each valid source grid cell using available monthly EWH values. Duplicate observations in one month are averaged; **missing months are not interpolated**. It then computes regional elastic displacement rates using the global source trend field.

**Important distinction:** The current repository uses **regular-grid load cells and centre-point Green's-function evaluation**. It does **not** contain the subsequently discussed adaptive Voronoi / MPAS source mesh, conservative polygon-overlap remapping, Bayesian inversion, or automatic error-controlled mesh optimization.

## Requirements

- MATLAB with support for `datetime`, `readmatrix`, `writetable`, and `exportgraphics`.
- For `mass_load_hv.m`: the function `distance` (typically from MATLAB Mapping Toolbox) and an external `calArea` implementation. **`calArea.m` is not currently included in this repository.**
- A compatible `cf_PREM.txt` file (provided).
- For the CSR example: a **separate** CSR mascon MAT file and a polygon boundary text file. Neither is included in this repository.

The CSR example itself implements its own spherical cell-area calculation and distance/azimuth helper, rather than calling `calArea` or `distance`.

## Quick start A — loading at a set of locations

Prepare a regular latitude–longitude grid `LLZ` with:

- `LLZ.lon`: [nLat × nLon] grid of **longitudes in degrees**
- `LLZ.lat`: [nLat × nLon] grid of **latitudes in degrees**
- `LLZ.rg`: [nLat × nLon] EWH in **metres**
- `sites`: [nSite × 2] array with columns **[longitude, latitude] in degrees**

For example, assuming the required auxiliary functions are on the MATLAB path:

```matlab
addpath('/path/to/ElasticLoading');
cd('/path/to/ElasticLoading');  % green_hv_region reads cf_PREM.txt here

% Replace these placeholders with real, correctly oriented EWH data:
% [LON,LAT] = meshgrid(lon_vector,lat_vector);
% LLZ = struct('lon',LON,'lat',LAT,'rg',ewh_m);
sites = [114.2 22.3; 100.0 35.0];

UEN_mm = green_hv_region(LLZ, sites);
% Rows correspond to sites, columns: [Up, East, North] (mm)
```

You can also evaluate an individual point with:

```matlab
T = load('cf_PREM.txt');
Gfun.theta = T(:,1);
Gfun.ur2   = T(:,4);
Gfun.v2    = T(:,5);

rho = 1000;  % kg/m^3
UEN_mm = mass_load_hv(Gfun, LLZ, [114.2 22.3], rho);
```

**Numerical warning:** The standalone `mass_load_hv` function evaluates a centre-point approximation and divides by angular separation; coincident source/receiver points require special near-field treatment. The current implementation changes non-finite contributions to zero, which is **not** a validated substitute for near-field cell integration.

## Quick start B — CSR mascon trends to regional displacement rates

Open [`P08_CSR_mascon_elastic_seafloor_PREM.m`](P08_CSR_mascon_elastic_seafloor_PREM.m) and edit the **SETTINGS** section, particularly:

```matlab
root      = '/path/to/your/research_project';
srcFile   = fullfile(root,'data','processed','GRACE','GRACE_mascon_4_20260707.mat');
boundFile = fullfile(root,'config','regions','NorthPacific_160E_130W_35N.txt');
greenFile = fullfile(root,'config','greens','cf_PREM.txt');

sourceEWHUnit = 'cm';  % MUST verify against original CSR data metadata
targetStepDeg = 2;     % output sampling resolution, degrees
nearAngleDeg  = 0.25;  % angular floor used in the Green's-function lookup
```

Expected CSR structure:

```text
mascon_csr.lon   longitudes (vector or 2-D grid)
mascon_csr.lat   latitudes  (vector or 2-D grid)
mascon_csr.rg    EWH, [nLat × nLon × nTime]
mascon_csr.tt    time: datetime or Nx3 [decimal year, year, month]
```

The boundary file must contain **[longitude, latitude]** columns, with optional NaN-separated polygon rings. The example uses a common calendar from **2005-01 through 2024-12**, with `minTrendMonths = 100` by default.

Run from MATLAB:

```matlab
OUT = P08_CSR_mascon_elastic_seafloor_PREM();
```

The configured output folders contain:

- `P08_CSR_elastic_PREM_results.mat`: calculated fields, settings, and summary values in `OUT`.
- `P08_CSR_elastic_displacement_sites.csv`: longitude, latitude, and Up/East/North rates (mm/yr).
- `CSR_elastic_up_trend_xyz.txt`: GMT-compatible point values.
- `CSR_EWH_trend_global_xyz.txt`: global source EWH trends.
- `P08_CSR_elastic_up_map.png` / `.pdf`: regional vertical-displacement-rate map (the current repository version uses a scatter plot).

The script also prints the area-weighted regional mean displacement rates. The output point spacing is an **evaluation-grid setting**, not the native GRACE spatial resolution.

## Limitations and recommended checks

1. **Input units:** Verify whether the particular CSR EWH field is in m, cm, or mm. An incorrect conversion directly scales the inferred load and displacement.
2. **Near-field integration:** The CSR example applies an angular floor (`nearAngleDeg`) to regularize the point-centre response. Results may be sensitive to this choice; it is not an exact finite-area load integral.
3. **Spatial resolution:** GRACE mascon gridded EWH values are not equivalent to fully independent native-resolution source masses. Interpret fine-scale signals carefully.
4. **Corrections:** The workflow does not automatically perform GIA correction, leakage correction, earthquake-step treatment, or ocean/land signal separation. Verify the input product's preprocessing before interpreting the result.
5. **Trend fitting:** The provided example uses an unweighted linear trend per grid cell and does not estimate serial-correlation-corrected trend uncertainties.
6. **Validation:** Compare source-area integration, numerical convergence with target spacing and near-field treatment, and independent reference calculations before publication.

## Citation and attribution

If this toolkit supports a publication, please cite the repository and the underlying loading-theory and Earth-model sources appropriate to your application. The present repository does **not** provide a versioned DOI or a formal software citation.

The load-convolution approach follows the classical surface-loading Green's-function framework associated with Farrell (1972); the provided coefficient table is labelled PREM. Users should verify the provenance and normalization of the tabulated Green's functions when preparing a publication.

## License

Distributed under the [MIT License](LICENSE). See the LICENSE file for the repository's copyright notice.

## Contact

**Chistrong Wen**  
The Hong Kong Polytechnic University  
Email: chistrong.wen@polyu.edu.hk

Contributions, reproducible test cases, and issue reports are welcome.
