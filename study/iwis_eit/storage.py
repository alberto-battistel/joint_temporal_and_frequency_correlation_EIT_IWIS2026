"""File paths, native lookup caches, checkpoints and metadata for the three scripts."""
from dataclasses import dataclass
import hashlib
import importlib.metadata
import json
from pathlib import Path
import sys

import numpy as np

from .monte_carlo import digest, load_condition, require_compatible, save_condition
from .numerics import endpoint_image


@dataclass
class DeterministicInput:
    arrays: dict
    pixel_dofs: np.ndarray
    metadata: dict
    directory: Path


def save_deterministic(directory, arrays, meshes, mapping, *, distance, mode,
                       penalty, forward_inverse, selected_frames):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    mesh_hashes = {}
    for name, mesh in meshes.items():
        path = directory / (name + '.vol')
        mesh.ngmesh.Save(str(path))
        mesh_hashes[name] = digest(path)
    coarse = meshes['coarse']
    arrays = dict(arrays,
                  coarse_points=np.array([v.point[:2] for v in coarse.vertices]),
                  coarse_cells=np.array([[v.nr for v in el.vertices] for el in coarse.Elements()]),
                  element_to_dof=np.array([mapping.space.GetDofNrs(el)[0] for el in coarse.Elements()]))
    np.savez_compressed(directory / 'deterministic.npz', **arrays)
    metadata = dict(
        schema_version=1, distance=distance, mode=mode, added_noise=0, lambda_value=penalty,
        lambda_convention='||M*x-y||^2 + lambda*||x||^2',
        axes={key: list(value.shape) for key, value in arrays.items()},
        voltage_axes=['injection-major electrode channel', 'time', 'frequency'],
        basis_rows=['static', 'frequency', 'cos', 'sin', 'frequency*cos', 'frequency*sin'],
        selected_frames=[list(frame) for frame in selected_frames],
        potential_reference='current electrodes mean', channel_count=256,
        include_current_electrodes=True, forward_inverse=forward_inverse,
        adjoint_inverse='sparsecholesky', mesh_elements={name: mesh.ne for name, mesh in meshes.items()},
        mesh_sha256=mesh_hashes, python=sys.version,
        units={'geometry': 'unspecified; metres if SI interpretation', 'conductivity': 'S/m',
               'current': '+/-1; A/m under SI unit-depth interpretation', 'voltage': 'V',
               'contact_impedance': 'ohm m^2 under SI interpretation',
               'thickness': 'unspecified; no conversion applied'},
        versions={name: importlib.metadata.version(name)
                  for name in ['ngsolve', 'netgen-mesher', 'numpy', 'scipy', 'matplotlib']})
    (directory / 'deterministic.json').write_text(json.dumps(metadata, indent=2) + '\n')

    # Cache the native P0 parameter index at each endpoint pixel, including its outside mask.
    lookup_field = mapping.to_field(np.arange(mapping.ndof) + 1)
    lookup_image = endpoint_image(lookup_field, coarse)
    pixel_dofs = np.full((64, 64), -1, dtype=int)
    valid = np.isfinite(lookup_image)
    pixel_dofs[valid] = np.rint(lookup_image[valid]).astype(int) - 1
    np.save(directory / 'pixel_dofs.npy', pixel_dofs)


def load_deterministic(project_root, position, mode):
    directory = Path(project_root) / 'results' / f'deterministic_d{position:.1f}_{mode}_v1'
    if not (directory / 'pixel_dofs.npy').exists():
        raise FileNotFoundError(f'Missing {directory}/pixel_dofs.npy. Run script 01 with this DISTANCE and MODE first.')
    with np.load(directory / 'deterministic.npz', allow_pickle=False) as saved:
        arrays = {key: saved[key] for key in saved.files}
    metadata = json.loads((directory / 'deterministic.json').read_text())
    if metadata['mode'] != mode or metadata['distance'] != position:
        raise ValueError('Wrong deterministic input')
    if arrays['clean_voltage'].shape != (256, 20, 3) or arrays['basis'].shape != (6, 60):
        raise ValueError('Unexpected deterministic axes')
    lookup = np.load(directory / 'pixel_dofs.npy', allow_pickle=False)
    if lookup.shape != (64, 64) or lookup.min() < -1 or lookup.max() >= arrays['jacobian'].shape[1]:
        raise ValueError('Invalid pixel-to-P0-DOF map')
    return DeterministicInput(arrays, lookup, metadata, directory)


def condition_directory(project_root, position, noise, mode, penalty):
    return (Path(project_root) / 'results/monte_carlo' / mode /
            f'd{position:.1f}_noise{noise:g}_lambda{penalty:g}_v2')


def save_monte_carlo(project_root, baseline, arrays, *, position, noise, mode,
                     penalty, seed, trial_count, representatives):
    data = baseline.arrays
    metadata = dict(
        schema_version=2, algorithm='paired-tikhonov-v2', mode=mode,
        position=position, noise_level=noise,
        noise_formula='eta * std(clean_voltage.ravel(), ddof=0) * N(0,1)',
        noise_amplitude=float(noise * np.std(data['clean_voltage'].ravel())),
        trial_count=trial_count, parameter_count=data['jacobian'].shape[1], seed=seed,
        representative_trials=list(representatives),
        representative_policy='fixed zero-based first three trials',
        lambda_value=penalty, lambda_convention='||M*x-y||^2 + lambda*||x||^2',
        rng_policy='per condition: SeedSequence(seed).spawn(n); child.generate_state(1,uint64)[0]; default_rng(child uint64); standard_normal(X.shape) once per trial; same tensor both methods',
        metric_policy='absolute Pearson on jointly finite endpoint pixels; undefined=NaN; median and linear nanquantiles for each method',
        sampling=dict(times_s=data['times_s'].tolist(), frequencies_hz=data['frequencies_hz'].tolist()),
        protocol='16 adjacent +/-1 injections (A/m under SI unit-depth interpretation); 16 single-ended electrodes per injection; includes driven electrodes',
        reference='CURRENT_ELECTRODES_MEAN', basis_sha256=hashlib.sha256(data['basis'].tobytes()).hexdigest(),
        axes=dict(metrics=['trial', 'circle1/circle2/circle3'],
                  images=['representative trial', 'component/inclusion', 'y', 'x'],
                  parameters=['representative trial', 'P0 DOF', 'component/inclusion']),
        deterministic_configuration=baseline.metadata,
        baseline_npz_sha256=digest(baseline.directory / 'deterministic.npz'),
        sampling_cache_sha256=digest(baseline.directory / 'pixel_dofs.npy'),
        deterministic_input=str(baseline.directory.relative_to(project_root)),
        python=sys.version, numpy_version=np.__version__)
    directory = condition_directory(project_root, position, noise, mode, penalty)
    _, _, summaries = save_condition(directory, arrays, metadata)
    return directory, summaries


def load_plot_conditions(project_root, positions, noises, mode, penalty):
    project_root = Path(project_root)
    conditions, input_hashes = [], {}
    for position in positions:
        for noise in noises:
            directory = condition_directory(project_root, position, noise, mode, penalty)
            if not (directory / 'trials.npz').exists():
                raise FileNotFoundError(f'Missing {directory}/trials.npz. Run scripts 01 and 02 for this configuration first.')
            data, metadata, _ = load_condition(directory)
            if (metadata['position'] != position or metadata['noise_level'] != noise or
                    metadata['lambda_value'] != penalty or metadata['mode'] != mode):
                raise ValueError(f'Selection does not match metadata: {directory}')
            baseline = project_root / metadata['deterministic_input']
            if digest(baseline / 'deterministic.npz') != metadata['baseline_npz_sha256']:
                raise ValueError('Deterministic input identity changed; use matching result inputs')
            with np.load(baseline / 'deterministic.npz', allow_pickle=False) as saved:
                geometry = {key: saved[key] for key in ['inclusion_centers', 'inclusion_radii']}
            for path in [directory / 'trials.npz', directory / 'metadata.json',
                         directory / 'summary.json', baseline / 'deterministic.npz']:
                input_hashes[str(path.relative_to(project_root))] = digest(path)
            conditions.append((directory.name, data, metadata, geometry))
    require_compatible([condition[2] for condition in conditions])
    return conditions, input_hashes


def save_export_metadata(project_root, directory, input_hashes, exported, config, *,
                         mode, representative_trial, save_figures, show_figures):
    if any(digest(Path(project_root) / path) != identity for path, identity in input_hashes.items()):
        raise RuntimeError('A scientific input changed during plotting')
    metadata = dict(
        mode=mode, input_sha256=input_hashes, png_dpi=config['png_dpi'], figure_config=config,
        svg_text='editable text', save_figures=save_figures, show_figures=show_figures,
        svg_images='eitlab.plotting.plot_image; embedded pixel images with vector text and overlays',
        representative_trial=representative_trial,
        scale_policy='independent max-absolute panel normalization; common [-1,1] colormap; amplitude comparisons not implied',
        boxplot_policy='Single/Unmixed Q1/Q3, median and 1.5 IQR whiskers',
        combined_boxplot_grouping='circle, then selected noise order, then Single/Unmixed',
        colormap='viridis', outputs=[path.name for path in exported] + ['summary.csv', 'summary.json', 'summary.md'])
    (Path(directory) / 'export_metadata.json').write_text(json.dumps(metadata, indent=2) + '\n')
