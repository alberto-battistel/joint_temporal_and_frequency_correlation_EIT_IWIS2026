"""Figure-only JSON settings; uses Python's standard library."""
import json
import math
from pathlib import Path

from matplotlib.text import Text


def load_figure_config(path):
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f'Missing figure settings: {path}. Use the supplied figure_config.json.')
    config = json.loads(path.read_text())

    def positive(value, name):
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
            raise ValueError(f'{name} must be a positive number')

    positive(config['png_dpi'], 'png_dpi')
    limits = config.get('boxplot_y_limits', [0, 1.03])
    if isinstance(limits, dict) and any(key not in ('0.4', '0.6') for key in limits):
        raise ValueError('boxplot_y_limits distance keys must be "0.4" or "0.6"')
    ranges = limits.items() if isinstance(limits, dict) else [('all', limits)]
    for distance, bounds in ranges:
        if bounds is None:
            continue
        if (not isinstance(bounds, list) or len(bounds) != 2 or
                any(isinstance(v, bool) or not isinstance(v, (int, float)) or
                    not math.isfinite(v) for v in bounds) or bounds[0] >= bounds[1]):
            raise ValueError(f'boxplot_y_limits ({distance}) must be [minimum, maximum] with finite minimum < maximum, or null')
    for key in ['base', 'titles', 'labels', 'ticks', 'legend', 'annotations']:
        positive(config['fonts'][key], f'fonts.{key}')
    sizes = {key: config['figure_sizes'][key] for key in
             ['geometry', 'conductivity', 'voltages', 'boxplot', 'combined_boxplot']}
    sizes.update({key: config['reconstruction'][key] for key in ['horizontal_size', 'vertical_size']})
    for key, size in sizes.items():
        if not isinstance(size, list) or len(size) != 2:
            raise ValueError(f'{key} must be [width, height] in centimetres')
        for value in size:
            positive(value, key)
    if config['reconstruction']['layout'] not in ('horizontal', 'vertical'):
        raise ValueError('reconstruction.layout must be horizontal or vertical')
    reconstruction = config['reconstruction']
    if not isinstance(reconstruction.get('compact', False), bool):
        raise ValueError('reconstruction.compact must be true or false')
    positive(reconstruction.get('panel_gap_cm', 0.15), 'reconstruction.panel_gap_cm')
    return config


def circle_label(name):
    """Readable circle labels without changing saved or geometry identifiers."""
    name = str(name)
    return 'circle ' + name[6:] if name.startswith('circle') and name[6:].isdigit() else name


def boxplot_y_limits(config, position):
    """Select limits for this distance; a missing distance uses automatic limits."""
    limits = config.get('boxplot_y_limits', [0, 1.03])
    return limits.get(f'{float(position):g}') if isinstance(limits, dict) else limits


def figure_size_inches(size_cm):
    """Convert configured [width, height] in centimetres for Matplotlib."""
    return tuple(value / 2.54 for value in size_cm)


def configure_figure(fig, config):
    """Apply font sizes to this figure only; leave global Matplotlib defaults alone."""
    fonts = config['fonts']
    for text in fig.findobj(match=Text):
        text.set_fontsize(fonts['base'])
    for ax in fig.axes:
        for title in (ax.title, ax._left_title, ax._right_title):
            title.set_fontsize(fonts['titles'])
        ax.xaxis.label.set_fontsize(fonts['labels'])
        ax.yaxis.label.set_fontsize(fonts['labels'])
        ax.tick_params(axis='both', which='both', labelsize=fonts['ticks'])
        ax.xaxis.get_offset_text().set_fontsize(fonts['ticks'])
        ax.yaxis.get_offset_text().set_fontsize(fonts['ticks'])
        for text in ax.texts:
            text.set_fontsize(fonts['annotations'])
        legend = ax.get_legend()
        if legend is not None:
            for text in [*legend.get_texts(), legend.get_title()]:
                text.set_fontsize(fonts['legend'])
    for text in fig.texts:
        text.set_fontsize(fonts['annotations'])
    if fig._suptitle is not None:
        fig._suptitle.set_fontsize(fonts['titles'])
