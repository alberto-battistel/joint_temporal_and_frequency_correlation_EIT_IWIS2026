"""Monte Carlo: visible noise, Single inversion, unmixing and image comparison."""

# %% Imports and editable parameters
from pathlib import Path
import numpy as np
from eitlab import TikhonovSolver
from iwis_eit.monte_carlo import cache_inverse_operator, correlation, run_seeds, TrialArchive
from iwis_eit.numerics import cached_parameter_images
from iwis_eit.storage import load_deterministic, save_monte_carlo

PROJECT_ROOT = Path(__file__).resolve().parents[1] if '__file__' in globals() else Path.cwd().parent
MODE = 'full'  # 'full': both distances, three noises, 100 trials per condition.
POSITIONS = (0.4,) if MODE == 'demo' else (0.4, 0.6)
NOISE_LEVELS = (0.1,) if MODE == 'demo' else (0.1, 0.01, 0.001)
SEED = 20260922
PENALTY = 1e-4  # Minimize ||M x - y||² + λ ||x||².
TRIAL_COUNT = {'demo': 3, 'full': 100}[MODE]
REPRESENTATIVE_TRIALS = (0, 1, 2)  # Fixed first three; never selected by appearance.


# %% One trial: the complete scientific reconstruction, one function deep
def reconstruct_trial(
    voltages,
    reference_voltage,
    inverse_operator,
    basis_pinv,
    pixel_dofs,
    reference_images,
    selected_frames,
):
    """Single differences and absolute-voltage unmixing use the SAME noisy tensor."""
    # Single: pick one time/frequency frame per circle and subtract the reference.
    single_voltage = np.column_stack(
        [voltages[:, t, f] - reference_voltage for t, f in selected_frames]
    )
    single_parameters = inverse_operator @ single_voltage

    # Unmixed: Y ≈ B D -> B = Y D⁺; invert all six voltage coefficients.
    fitted_voltage = voltages.reshape(voltages.shape[0], -1) @ basis_pinv
    unmixed_parameters = inverse_operator @ fitted_voltage

    # Apply the native pixel lookup saved by 01; this does not solve or resample FEM.
    single_images = cached_parameter_images(single_parameters, pixel_dofs)
    component_images = cached_parameter_images(unmixed_parameters, pixel_dofs)
    unmixed_images = np.stack(
        [
            component_images[1],
            np.hypot(component_images[2], component_images[3]),
            np.hypot(component_images[4], component_images[5]),
        ]
    )
    single_correlations = np.array(
        [correlation(reference, image) for reference, image in zip(reference_images, single_images)]
    )
    unmixed_correlations = np.array(
        [
            correlation(reference, image)
            for reference, image in zip(reference_images, unmixed_images)
        ]
    )
    return dict(
        single_parameters=single_parameters,
        unmixed_parameters=unmixed_parameters,
        single_images=single_images,
        component_images=component_images,
        unmixed_images=unmixed_images,
        single_correlations=single_correlations,
        unmixed_correlations=unmixed_correlations,
    )


# %% Load inputs; create an EITLab inverse and unmixing operator once per distance
prepared = {}
for position in POSITIONS:
    baseline = load_deterministic(PROJECT_ROOT, position, MODE)
    data = baseline.arrays
    solver = TikhonovSolver(PENALTY, data['jacobian'])
    # Cache the linear solver's action once; each trial then uses matrix products.
    inverse_operator = cache_inverse_operator(solver, data['jacobian'].shape[0])
    basis_pinv = np.linalg.pinv(data['basis'])
    prepared[position] = (baseline, inverse_operator, basis_pinv)

# %% Noise and paired trials: all scientific operations are visible here or above
for position in POSITIONS:
    baseline, inverse_operator, basis_pinv = prepared[position]
    data = baseline.arrays
    clean_voltage = data['clean_voltage']
    selected_frames = baseline.metadata['selected_frames']
    for noise_level in NOISE_LEVELS:
        noise_amplitude = noise_level * np.std(clean_voltage.ravel())
        seeds = run_seeds(SEED, TRIAL_COUNT)
        archive = TrialArchive(TRIAL_COUNT, data['inclusion_names'], seeds, REPRESENTATIVE_TRIALS)
        for trial, child_seed in enumerate(seeds):
            rng = np.random.default_rng(child_seed)
            noise = noise_amplitude * rng.standard_normal(clean_voltage.shape)
            noisy_voltage = clean_voltage + noise
            result = reconstruct_trial(
                noisy_voltage,
                data['reference_voltage'],
                inverse_operator,
                basis_pinv,
                baseline.pixel_dofs,
                data['reference_images'],
                selected_frames,
            )
            archive.record(trial, noise, result)

        # Storage/metadata bookkeeping is deliberately outside the scientific flow.
        output, summaries = save_monte_carlo(
            PROJECT_ROOT,
            baseline,
            archive.arrays(),
            position=position,
            noise=noise_level,
            mode=MODE,
            penalty=PENALTY,
            seed=SEED,
            trial_count=TRIAL_COUNT,
            representatives=REPRESENTATIVE_TRIALS,
        )
        print(
            'Saved',
            output.name,
            TRIAL_COUNT,
            'trials; median |r|:',
            'Single',
            summaries['single']['median'],
            'Unmixed',
            summaries['unmixed']['median'],
            flush=True,
        )
