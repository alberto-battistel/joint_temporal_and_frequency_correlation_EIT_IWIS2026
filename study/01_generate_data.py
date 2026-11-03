"""Generate deterministic study inputs; use 02 for trials and 03 for figures. Run sequential cells from study/ or set PROJECT_ROOT.

Requires an installed compatible EITLab. Set MODE and DISTANCE below.
"""

# %% Imports and editable configuration
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from eitlab import (
    ElectrodeDisk,
    ElectricalModel,
    ElectricalLinearization,
    JacobianAssembler,
    MeasurementPlan,
    MeasurementRunner,
    ProtocolFactory,
    SingleEnded,
    PotentialReference,
    ElementwiseConductivity,
    TikhonovSolver,
)
from iwis_eit.numerics import endpoint_image
from iwis_eit.monte_carlo import correlation
from iwis_eit.storage import save_deterministic
from iwis_eit.publication import export_figure
from iwis_eit.figure_config import (
    load_figure_config,
    configure_figure,
    figure_size_inches,
    circle_label,
)

PROJECT_ROOT = Path(__file__).resolve().parents[1] if '__file__' in globals() else Path.cwd().parent
FIGURE_CONFIG = load_figure_config(PROJECT_ROOT / 'figure_config.json')
DISTANCE = 0.6
MODE = 'full'  # 'full': maxh=0.1/0.05; 'demo': maxh=0.22/0.16
OUTPUT_DIRECTORY = PROJECT_ROOT / 'results' / f'deterministic_d{DISTANCE:.1f}_{MODE}_v1'
PENALTY = 1e-4  # Minimize ||M x - y||² + λ ||x||².
SELECTED_FRAMES = ((0, 0), (2, 2), (4, 2))  # One (time, frequency) index per circle.
FORWARD_INVERSE = 'sparsecholesky'  # EITLab forwards this to its FEM backend.
SAVE_FIGURES = True  # Write SVG and PNG files.
SHOW_FIGURES = True  # Display figures independently of saving.

# %% Geometry: disk, 16 electrodes, three inclusions
centers = np.array(
    [
        (DISTANCE * np.cos((i + 1) * 2 * np.pi / 3), DISTANCE * np.sin((i + 1) * 2 * np.pi / 3))
        for i in range(3)
    ]
)
names = ['circle1', 'circle2', 'circle3']
coarse_geometry = ElectrodeDisk(radius=1, n_electrodes=16)
homogeneous_geometry = ElectrodeDisk(radius=1, n_electrodes=16)
inclusion_geometry = ElectrodeDisk(radius=1, n_electrodes=16)
for name, center in zip(names, centers):
    inclusion_geometry.add_disk_inclusion(name=name, p=tuple(center), r=0.3)
if MODE not in ('demo', 'full'):
    raise ValueError('MODE must be demo or full')

# %% Meshing: coarse inverse mesh and independent forward meshes
coarse_mesh = coarse_geometry.make_mesh(maxh=0.1 if MODE == 'full' else 0.22)
homogeneous_mesh = homogeneous_geometry.make_mesh(maxh=0.05 if MODE == 'full' else 0.16)
inclusion_mesh = inclusion_geometry.make_mesh(maxh=0.05 if MODE == 'full' else 0.16)

# %% Measurement protocol: adjacent current injections, single-ended voltages
electrodes = tuple(name for name in coarse_mesh.GetBoundaries() if name.startswith('el_'))
protocol = ProtocolFactory.adjacent(electrodes, amplitude=1)
plan = MeasurementPlan(protocol, SingleEnded)

# %% EITLab forward models and conductivity parameters
model_options = dict(
    electrode_names=electrodes,
    numerical_ground_el=electrodes[0],
    z_contact=1e-3,
    potential_reference=PotentialReference.CURRENT_ELECTRODES_MEAN,
    order=2,
    inverse=FORWARD_INVERSE,
    complex_mode=False,
)
coarse_model = ElectricalModel(coarse_mesh, sigma=1, **model_options)
homogeneous_model = ElectricalModel(homogeneous_mesh, sigma=1, **model_options)
inclusion_model = ElectricalModel(inclusion_mesh, sigma=1, **model_options)
mapping = ElementwiseConductivity(coarse_mesh)

# %% Reference forward solutions and adjoint Jacobian M
coarse_result = MeasurementRunner(coarse_model, plan).run()
linearization = ElectricalLinearization(coarse_model, mapping, adjoint_inverse='sparsecholesky')
M = JacobianAssembler(linearization, coarse_result).assemble_total_jacobian()
X0 = MeasurementRunner(homogeneous_model, plan).run(store_solutions=False).vector

# %% Time/frequency conductivity laws and absolute voltages, frequency fastest
# Full documented sampling is retained even in the coarse demonstration.
times = np.arange(20) / 20 * 2
frequencies = np.array([1e3, 1e4, 1e5])
f_normalized = (frequencies - frequencies.mean()) / frequencies.std()
# Circle 1: linear spectral law with specified endpoint conductivities (S/m).
sigma1_min, sigma1_max = 0.1, 0.6
frequency_fraction = (f_normalized - f_normalized.min()) / (f_normalized.max() - f_normalized.min())
sigma1_frequency = sigma1_min + (sigma1_max - sigma1_min) * frequency_fraction
# Conductivities in S/m: spectral / temporal / combined components.
sigmas = np.stack(
    [
        np.broadcast_to(sigma1_frequency, (20, 3)),
        np.broadcast_to((0.5 + 0.2 * np.cos(4 * np.pi * times + 1))[:, None], (20, 3)),
        0.4 + 0.2 * np.cos(4 * np.pi * times + 1 + np.pi)[:, None] * f_normalized,
    ]
)
X = np.empty((X0.size, len(times), len(frequencies)))
for t in range(len(times)):
    for f in range(len(frequencies)):
        inclusion_model.set_sigma(
            dict(body=1, **{name: sigmas[i, t, f] for i, name in enumerate(names)})
        )
        X[:, t, f] = MeasurementRunner(inclusion_model, plan).run(store_solutions=False).vector

# %% Reference spatial images used for the absolute Pearson correlations
# These are fixed shape templates (0.1 S/m inside one circle, 1 S/m elsewhere).
reference_images = []
for name in names:
    inclusion_model.set_sigma(
        {material: (0.1 if material == name else 1.0) for material in inclusion_mesh.GetMaterials()}
    )
    reference_images.append(endpoint_image(inclusion_model.sigma_cf, inclusion_mesh))
reference_images = np.stack(reference_images)

# %% Unmixing basis D: constant, frequency, temporal harmonics and their products
# Frequency is the fastest axis in X.reshape(channel_count, time_count * frequency_count).
cos_time = np.cos(4 * np.pi * times)
sin_time = np.sin(4 * np.pi * times)
f = np.tile(f_normalized, len(times))
c = np.repeat(cos_time, len(frequencies))
s = np.repeat(sin_time, len(frequencies))
D = np.vstack([np.ones(f.size), f, c, s, c * f, s * f])

# %% Inversion: EITLab Tikhonov solver and the three Single voltage differences
solver = TikhonovSolver(PENALTY, M)
single_voltage = np.column_stack([X[:, t, f] - X0 for t, f in SELECTED_FRAMES])
single_parameters = np.column_stack([solver.solve(column) for column in single_voltage.T])

# %% Unmix absolute voltages first, then invert each of the six coefficient columns
# Y ≈ B D -> B = Y D⁺. Do not subtract X0 from the absolute-voltage tensor here.
fitted_voltage = X.reshape(X.shape[0], -1) @ np.linalg.pinv(D)
unmixed_parameters = np.column_stack([solver.solve(column) for column in fitted_voltage.T])

# %% Images: EITLab parameter fields sampled on the fixed 64 x 64 endpoint grid
single_images = []
for parameters in single_parameters.T:
    field = mapping.to_field(parameters, positive=False)
    single_images.append(endpoint_image(field, coarse_mesh))
single_images = np.stack(single_images)

component_images = []
for parameters in unmixed_parameters.T:
    field = mapping.to_field(parameters, positive=False)
    component_images.append(endpoint_image(field, coarse_mesh))
component_images = np.stack(component_images)
# Compare the spectral coefficient and the two phase-independent harmonic norms.
unmixed_images = np.stack(
    [
        component_images[1],
        np.hypot(component_images[2], component_images[3]),
        np.hypot(component_images[4], component_images[5]),
    ]
)
single_correlations = np.array([correlation(r, v) for r, v in zip(reference_images, single_images)])
unmixed_correlations = np.array(
    [correlation(r, v) for r, v in zip(reference_images, unmixed_images)]
)

# %% Save study arrays; native checkpoints, pixel caching and metadata stay in the helper
study_arrays = dict(
    jacobian=M,
    reference_voltage=X0,
    clean_voltage=X,
    jacobian_reference_voltage=coarse_result.vector,
    basis=D,
    times_s=times,
    frequencies_hz=frequencies,
    frequency_normalized=f_normalized,
    inclusion_sigma=sigmas,
    inclusion_centers=centers,
    inclusion_radii=np.full(3, 0.3),
    inclusion_names=np.array(names),
    electrode_names=np.array(electrodes),
    reference_images=reference_images,
    single_parameters=single_parameters,
    unmixed_parameters=unmixed_parameters,
    fitted_voltage=fitted_voltage,
    single_images=single_images,
    component_images=component_images,
    unmixed_images=unmixed_images,
    single_correlations=single_correlations,
    unmixed_correlations=unmixed_correlations,
)
save_deterministic(
    OUTPUT_DIRECTORY,
    study_arrays,
    meshes=dict(coarse=coarse_mesh, homogeneous=homogeneous_mesh, inclusion=inclusion_mesh),
    mapping=mapping,
    distance=DISTANCE,
    mode=MODE,
    penalty=PENALTY,
    forward_inverse=FORWARD_INVERSE,
    selected_frames=SELECTED_FRAMES,
)
print(
    'Saved',
    OUTPUT_DIRECTORY,
    'Single |r|:',
    single_correlations,
    'Unmixed |r|:',
    unmixed_correlations,
)

# %% Circle schematic: one Matplotlib cycle color per circle
circle_colors = plt.rcParams['axes.prop_cycle'].by_key()['color'][: len(names)]
geometry_fig, geometry_ax = plt.subplots(
    figsize=figure_size_inches(FIGURE_CONFIG['figure_sizes']['geometry']), layout='constrained'
)
domain_radius = coarse_geometry.radius
geometry_ax.add_patch(plt.Circle((0, 0), domain_radius, fill=False, color='black'))
# Use the electrode locations and arc widths from the EITLab geometry.
for index, electrode_name in enumerate(coarse_geometry.electrode_names):
    start_angle = coarse_geometry.angular_shift + index * coarse_geometry.sector_angle
    angles = np.deg2rad(np.linspace(start_angle, start_angle + coarse_geometry.electrode_angle, 32))
    geometry_ax.plot(
        domain_radius * np.cos(angles),
        domain_radius * np.sin(angles),
        color='0.25',
        linewidth=4,
        solid_capstyle='butt',
        zorder=3,
    )
    middle_angle = np.deg2rad(start_angle + coarse_geometry.electrode_angle / 2)
    geometry_ax.text(
        1.14 * domain_radius * np.cos(middle_angle),
        1.14 * domain_radius * np.sin(middle_angle),
        electrode_name,
        ha='center',
        va='center',
        fontsize=7.5,
        color='0.25',
    )
for name, center, color in zip(names, centers, circle_colors):
    geometry_ax.add_patch(plt.Circle(center, 0.3, fill=False, color=color, linewidth=2))
    geometry_ax.text(*center, circle_label(name), ha='center', va='center', color=color)
geometry_ax.set(
    xlim=(-1.28 * domain_radius, 1.28 * domain_radius),
    ylim=(-1.28 * domain_radius, 1.28 * domain_radius),
    title=f'Circle schematic | d={DISTANCE:g} m | {MODE}',
)
geometry_ax.set_aspect('equal')
geometry_ax.set_axis_off()

# %% Conductivity laws: one subplot per circle, with matching schematic colors
conductivity_fig, conductivity_axes = plt.subplots(
    3,
    1,
    figsize=figure_size_inches(FIGURE_CONFIG['figure_sizes']['conductivity']),
    sharex=True,
    sharey=True,
    layout='constrained',
)
frequency_styles = ['-', '--', ':']
for i, (name, ax, color) in enumerate(zip(names, conductivity_axes, circle_colors)):
    for k, (frequency, linestyle) in enumerate(zip(frequencies, frequency_styles)):
        ax.plot(
            times, sigmas[i, :, k], color=color, linestyle=linestyle, label=f'{frequency/1e3:g} kHz'
        )
    ax.set_title(circle_label(name), loc='left', color=color)
    ax.set_ylabel('σ (S/m)')
    ax.legend(title='Frequency', fontsize=8, title_fontsize=8, loc='lower right')
    ax.grid(alpha=0.2)
conductivity_axes[-1].set_xlabel('Time (s)')
conductivity_fig.suptitle(f'σ | d={DISTANCE:g} m | {MODE}')

# %% Referenced voltages in their own figure, using Matplotlib's color cycle
voltage_fig, voltage_ax = plt.subplots(
    figsize=figure_size_inches(FIGURE_CONFIG['figure_sizes']['voltages']), layout='constrained'
)
voltage_ax.plot(X0, label='Homogeneous reference')
voltage_ax.plot(X[:, 0, 0], label='Target: t=0 s, f=1 kHz')
voltage_ax.set(
    xlabel='Injection-major channel',
    ylabel='Voltage (V)',
    title=f'Current-electrode-mean reference | d={DISTANCE:g} m | {MODE}',
)
voltage_ax.legend()
voltage_ax.grid(alpha=0.2)

# %% Export/display: controlled independently; file-writing details stay in the helper
for figure, filename in [
    (geometry_fig, 'geometry'),
    (conductivity_fig, 'conductivity'),
    (voltage_fig, 'voltages'),
]:
    if SAVE_FIGURES:
        export_figure(figure, OUTPUT_DIRECTORY / filename, FIGURE_CONFIG)
    else:
        # Apply the same configured fonts to display-only figures.
        configure_figure(figure, FIGURE_CONFIG)
    if SHOW_FIGURES:
        plt.show()
    plt.close(figure)
