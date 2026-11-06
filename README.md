# Joint Temporal and Frequency Correlation for Multi-Frequency Electrical Impedance Tomography

This repository contains the scripts used to generate the simulations, reconstructions, figures, and tables for the study presented at the **19th International Workshop on Impedance Spectroscopy, September 22–25, 2026**, and described in the manuscript below.

## Associated publication

**Joint Temporal and Frequency Correlation for Multi-Frequency Electrical Impedance Tomography**

- Authors: `<AUTHORS>`
- Publication status: `<PUBLICATION_STATUS>`
- Article DOI: pending

## Software versions and citation

The software versions used for the manuscript are:

| Software | Release tag | Version-specific Zenodo DOI |
|---|---|---|
| This study repository | `<STUDY_TAG>` | `<STUDY_VERSION_DOI>` |
| EITLab | `<EITLAB_TAG>` | `<EITLAB_VERSION_DOI>` |

Release tags and Zenodo DOIs will be added when the corresponding releases are published.
To reproduce the manuscript results, use these releases and follow the **Full experiment** instructions below. The default demo configuration is intended for a quick run.

If you use this work, please cite the relevant software releases listed above. Use the version-specific DOIs to identify the exact archived software used.

## Installation

Use the EITLab release tagged **<EITLAB_TAG>**.

Extract this package. From the `joint_temporal_and_frequency_correlation_EIT_IWIS2026` folder, run the following commands, replacing `<EITLAB_REPOSITORY_URL>` and `<EITLAB_TAG>` with the repository URL and release tag:

```bash
conda create --name iwis2026 python=3.12 pip
conda activate iwis2026
git clone --branch "<EITLAB_TAG>" --depth 1 "<EITLAB_REPOSITORY_URL>" ../eitlab-dev
python -m pip install -r requirements.txt
```

The requirements file installs that EITLab checkout with its simulation and plotting dependencies.
If you already have EITLab checked out at the specified tag, use its folder path in `requirements.txt` and skip cloning.

Activate the environment and open the `joint_temporal_and_frequency_correlation_EIT_IWIS2026` folder before running the scripts.

## Quick run

Run these commands in order:

```bash
python study/01_generate_data.py
python study/02_run_monte_carlo.py
python study/03_make_figures_tables.py
```

The default `MODE='demo'` uses coarse meshes, distance 0.4, noise level 0.1, and three Monte Carlo trials.
Settings are near the beginning of each script.

Each script starts independently and reads the files produced by the preceding script.
Summary tables are printed to the console and saved to files.

The scientific sequence is visible in these three scripts:

- **01:** geometry and meshing, EITLab forward models and Jacobian, conductivity laws, the six-row unmixing basis, Tikhonov inversion, parameter-to-image conversion, and geometry/conductivity/voltage plotting.
- **02:** EITLab inverse setup, explicit noise generation, and `reconstruct_trial()`.
This  function shows Single frame selection, voltage unmixing, inversion, image construction and correlation calculation.
- **03:** saved-result summaries and `plot_reconstruction_images()`.
This function calls EITLab's `plot_image()` and shows the panel arrangement, normalization, geometry outlines and shared colorbar.

Helpers in `study/iwis_eit/` handle file loading/saving, metadata, inverse-operator caching, trial storage, pixel lookup and common figure formatting.


## Full experiment

1. Set `MODE='full'` in all three scripts.
2. In `01_generate_data.py`, set `DISTANCE=0.4` and run it:

   ```bash
   python study/01_generate_data.py
   ```

3. Change `DISTANCE=0.6` and run script 01 again.
4. Run the reconstructions and generate the plots:

   ```bash
   python study/02_run_monte_carlo.py
   python study/03_make_figures_tables.py
   ```

The full experiment uses both distances, noise levels 0.1, 0.01, and 0.001, and 100 trials per condition.

## Outputs

Data is saved in `results/`. Script 01 saves its diagnostic figures alongside its data. Script 03 saves plots and summary tables in `outputs/`.

- Script 01: circle schematic, three conductivity subplots with matching circle colors, and a separate voltage figure.
- Script 03: Single/Unmixed correlation boxplots, one combined boxplot per distance, representative EIT images, and CSV, Markdown, and JSON tables.

Each full-mode combined boxplot has six boxes per circle: one Single/Unmixed pair per noise level.
Figures are exported as SVG and PNG.
EIT panels are normalized independently, so compare spatial patterns rather than amplitudes.

Reconstruction panels use EITLab's `eitlab.plotting.plot_image()` to display the saved pixel arrays. SVG files embed these pixel images while retaining editable text and vector outlines.

## Regenerate plots

Set `MODE` in script 03 to match your saved results, then run:

```bash
python study/03_make_figures_tables.py
```

This uses saved data without repeating simulations or reconstructions. Scripts 01 and 03 each have independent `SAVE_FIGURES` and `SHOW_FIGURES` flags near the top:

- `SAVE_FIGURES=True`, `SHOW_FIGURES=True`: save SVG/PNG files and display plots.
- `SAVE_FIGURES=False`, `SHOW_FIGURES=True`: display plots without saving plot files.
- `SAVE_FIGURES=True`, `SHOW_FIGURES=False`: save plot files without opening plot windows.
- Both `False`: neither save nor display plots.

These flags control plots only.
Script 01 still saves simulation data;
script 03 still saves summary tables and export metadata.
Disabling plot saving leaves existing plot files untouched.
Repeating a configuration with saving enabled overwrites its corresponding outputs.


## Figure configuration

Edit `figure_config.json` to set figure sizes, font sizes and PNG resolution. All sizes are `[width, height]` in centimetres (cm) and set the exported canvas size; fonts are in points. JSON uses Python's standard library, so no extra dependency is needed.

- `figure_sizes`: individual sizes for the geometry, conductivity, voltage, condition boxplot and combined boxplot figures.
- `fonts`: sizes for base text, titles, labels, ticks, legends and annotations.
- `png_dpi`: PNG export resolution.
- `boxplot_y_limits`: independent `[minimum, maximum]` vertical-axis limits for distances `"0.4"` and `"0.6"`. Each distance's limits apply to both its individual and combined correlation boxplots. Use `null` for a distance to give it automatic limits; omitted distances also use automatic limits. This changes only the displayed axis range, not the data or statistics.
- `reconstruction.layout`: `"horizontal"` gives **2 rows x 3 columns** (Single/Unmixed rows, circle columns), with the colorbar on the right; `"vertical"` gives **3 rows x 2 columns** (circle rows, Single/Unmixed columns), with a horizontal colorbar below the images.
- `reconstruction.horizontal_size` and `reconstruction.vertical_size`: figure size for each layout.
- `reconstruction.compact`: `true` uses shared headers and compact labels, and omits the condition heading and explanatory footer. Use this for publication columns. Set it to `false` for the detailed diagnostic presentation. Reconstruction images have no coordinate axis labels or tick numbers in either mode; the shared colorbar retains its scale.
- `reconstruction.panel_gap_cm`: space reserved between panels in compact mode. The default is `0.15` cm; smaller values leave more room for the images.

The compact reconstruction defaults are **8.5 x 14 cm** for a single-column vertical figure and **18.5 x 13 cm** for a double-column horizontal figure.
Width and height both matter: square image panels need enough height for all rows.
Increase the height if you want larger panels while keeping the same column width.
Font sizes remain in points; the figure is exported at the configured size without cropping the canvas.

For example, to use different vertical scales for the two distances:

```json
"boxplot_y_limits": {
  "0.4": [0.2, 0.8],
  "0.6": [0.4, 0.9]
}
```

For example, to put Single and Unmixed side by side for each circle:

```json
"reconstruction": {
  "layout": "vertical",
  "compact": true,
  "panel_gap_cm": 0.15,
  "horizontal_size": [18.5, 13],
  "vertical_size": [8.5, 12]
}
```

In compact figures, Single and Unmixed are the column headers in the vertical layout; circles are the column headers in the horizontal layout. Corner labels identify the other grouping. The shared colorbar shows each panel divided by its own maximum absolute coefficient, from -1 to 1; compare spatial patterns rather than amplitudes. Dotted circles mark the true inclusions. Single panels use `(t, f) = (0 s, 1 kHz)`, `(0.2 s, 100 kHz)` and `(0.4 s, 100 kHz)` for circles 1–3. Unmixed panels show the frequency coefficient, temporal harmonic norm and frequency-dependent harmonic norm, respectively. Images use a 64 x 64 grid over `[-1, 1]` m in each direction. The fixed representative trial and condition parameters remain in the figure metadata exported by script 03 and can be included in your article caption.

After changing these settings, rerun script 03 to regenerate the reconstruction images and boxplots from saved results. To update script 01's diagnostic figures without recomputing data in an existing session, reload `FIGURE_CONFIG` and rerun its plotting and export cells. Both scripts read the same configuration file. Existing figure files are replaced with the newly configured exports.
