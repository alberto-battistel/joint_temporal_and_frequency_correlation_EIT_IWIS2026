"""Plot saved results: EITLab image rendering is explicit below, with no new trials."""

# %% Imports and editable selection
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from eitlab.plotting import plot_image
from iwis_eit.figure_config import (
    load_figure_config,
    configure_figure,
    figure_size_inches,
)
from iwis_eit.publication import (
    table_rows,
    markdown_table,
    export_tables,
    boxplots,
    combined_boxplot,
    condition_title,
    export_figure,
)
from iwis_eit.storage import load_plot_conditions, save_export_metadata

PROJECT_ROOT = Path(__file__).resolve().parents[1] if '__file__' in globals() else Path.cwd().parent
FIGURE_CONFIG = load_figure_config(PROJECT_ROOT / 'figure_config.json')
MODE = 'full'
POSITIONS = (0.4,) if MODE == 'demo' else (0.4, 0.6)
NOISE_LEVELS = (0.1,) if MODE == 'demo' else (0.1, 0.01, 0.001)
PENALTY = 1e-4
REPRESENTATIVE_TRIAL = 0  # Fixed first trial; never chosen by appearance.
SAVE_FIGURES = True
SHOW_FIGURES = True
OUTPUT_DIRECTORY = PROJECT_ROOT / 'outputs' / MODE / f'lambda{PENALTY:g}_v2'


# %% Reconstruction figure: layout, EITLab plot_image and geometry outlines, one function deep
def plot_reconstruction_images(data, m, deterministic, trial=0, *, config):
    """Draw six saved images through EITLab; no resampling or inversion."""
    indices = data['representative_trial_indices']
    matches = np.flatnonzero(indices == trial)
    if not len(matches):
        raise ValueError(f'Trial {trial} not saved; available: {indices.tolist()}')
    idx = int(matches[0])
    layout = config['reconstruction']['layout']
    compact = config['reconstruction'].get('compact', False)
    shape = (2, 3) if layout == 'horizontal' else (3, 2)
    size = config['reconstruction'][layout + '_size']
    fig, axes = plt.subplots(
        *shape, figsize=figure_size_inches(size), layout='compressed' if compact else 'constrained'
    )
    if compact:
        pad = config['reconstruction'].get('panel_gap_cm', 0.15) / 2 / 2.54
        fig.get_layout_engine().set(w_pad=pad, h_pad=pad, wspace=0, hspace=0)
    panels = np.stack([data['single_images'][idx], data['unmixed_images'][idx]])
    times = m['sampling']['times_s']
    frequencies = m['sampling']['frequencies_hz']
    selected_frames = m['deterministic_configuration']['selected_frames']
    selections = [
        f'Single: t={times[t]:g} s, f={frequencies[f]/1e3:g} kHz' for t, f in selected_frames
    ]
    components = [
        'Unmixed: frequency coefficient',
        'Unmixed: sqrt(cos² + sin²)',
        'Unmixed: sqrt((f cos)² + (f sin)²)',
    ]
    for row in range(2):
        for col in range(3):
            ax = axes[row, col] if layout == 'horizontal' else axes[col, row]
            image = panels[row, col]
            scale = float(np.nanmax(abs(image)))
            normalized = image / scale if scale > 0 else image.copy()
            artist = plot_image(
                np.ma.masked_invalid(normalized),
                ax=ax,
                extent=(-1, 1, -1, 1),
                colorbar=False,
                cmap='viridis',
                vmin=-1,
                vmax=1,
                interpolation='nearest',
            )
            ax.add_patch(Circle((0, 0), 1, fill=False, color='gray', linewidth=0.7))
            # Mark the true inclusion corresponding to this panel.
            center = deterministic['inclusion_centers'][col]
            radius = deterministic['inclusion_radii'][col]
            ax.add_patch(
                Circle(center, radius, fill=False, color='red', linewidth=0.9, linestyle='--')
            )
            ax.set(xlim=(-1, 1), ylim=(-1, 1))
            if compact:
                method = ['Single', 'Unmixed'][row]
                if layout == 'vertical':
                    if col == 0:
                        ax.set_title(method)
                else:
                    ax.text(
                        0.01,
                        0.99,
                        method,
                        transform=ax.transAxes,
                        ha='left',
                        va='top',
                        bbox=dict(facecolor='white', edgecolor='none', pad=0.3),
                    )
            else:
                title = selections[col] if row == 0 else components[col]
                ax.set_title(title)
            ax.set_aspect('equal')
            ax.set_facecolor('#eeeeee')
            ax.set_axis_off()
    colorbar_options = dict(fraction=0.035, pad=0.025) if layout == 'vertical' else {}
    if compact:
        colorbar_options['ticks'] = [-1, 0, 1]
    colorbar = fig.colorbar(
        artist,
        ax=axes,
        shrink=0.95 if compact else 0.7,
        location='bottom' if layout == 'vertical' else 'right',
        label=(
            'Normalized coefficient'
            if compact
            else 'Coefficient / panel maximum absolute coefficient'
        ),
        **colorbar_options,
    )
    colorbar.solids.set_rasterized(False)
    if not compact:
        fig.suptitle(condition_title(m) + f'\nRepresentative trial {trial} (fixed before analysis)')
        note = 'Each panel normalized independently; compare spatial patterns, not amplitudes. Red dashed circle: true inclusion for this panel.\nSingle: signed conductivity change; Unmixed: frequency coefficient / harmonic norms. Endpoint sampling: 64 x 64.'
        fig.supxlabel(
            note + '\nGeometry in metres under the SI/unit-depth interpretation.', fontsize=8
        )
    return fig


# %% Load matching result files; paths, consistency checks and hashing stay in the helper
conditions, input_hashes = load_plot_conditions(
    PROJECT_ROOT, POSITIONS, NOISE_LEVELS, MODE, PENALTY
)
print(f'Loaded {len(conditions)} {MODE} conditions')

# %% Recompute numeric table from raw trials, display, and export
rows = [row for _, data, metadata, _ in conditions for row in table_rows(data, metadata)]
export_tables(rows, OUTPUT_DIRECTORY)
print(markdown_table(rows))

# %% Correlation boxplots and fixed representative spatial comparisons
exported = []
for name, data, metadata, display_geometry in conditions:
    correlation_figure = boxplots(data, metadata, config=FIGURE_CONFIG)
    reconstruction_figure = plot_reconstruction_images(
        data, metadata, display_geometry, REPRESENTATIVE_TRIAL, config=FIGURE_CONFIG
    )
    for figure, suffix in [
        (correlation_figure, 'correlations'),
        (reconstruction_figure, f'reconstruction_trial{REPRESENTATIVE_TRIAL:03d}'),
    ]:
        configure_figure(figure, FIGURE_CONFIG)
        if SAVE_FIGURES:
            exported += export_figure(
                figure, OUTPUT_DIRECTORY / (MODE + '_' + name + '_' + suffix), FIGURE_CONFIG
            )
        if SHOW_FIGURES:
            plt.show()
        plt.close(figure)

# %% One additional combined correlation boxplot per distance
for position in POSITIONS:
    selected = [
        (data, metadata) for _, data, metadata, _ in conditions if metadata['position'] == position
    ]
    combined_figure = combined_boxplot(selected, config=FIGURE_CONFIG)
    stem = OUTPUT_DIRECTORY / f'{MODE}_d{position:.1f}_lambda{PENALTY:g}_combined_correlations'
    configure_figure(combined_figure, FIGURE_CONFIG)
    if SAVE_FIGURES:
        exported += export_figure(combined_figure, stem, FIGURE_CONFIG)
    if SHOW_FIGURES:
        plt.show()
    plt.close(combined_figure)

# %% Save export choices; storage and source-file identity checks stay in the helper
save_export_metadata(
    PROJECT_ROOT,
    OUTPUT_DIRECTORY,
    input_hashes,
    exported,
    FIGURE_CONFIG,
    mode=MODE,
    representative_trial=REPRESENTATIVE_TRIAL,
    save_figures=SAVE_FIGURES,
    show_figures=SHOW_FIGURES,
)
print('Exported to', OUTPUT_DIRECTORY)
for path in exported:
    print(path.name)
