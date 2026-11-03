"""Inverse-operator cache, trial storage, metrics and result-file checks."""
import hashlib
import json
from pathlib import Path
import warnings
import numpy as np

NAMES = np.array(['circle1','circle2','circle3'])

def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def run_seeds(seed, count):
    return np.array([s.generate_state(1,dtype=np.uint64)[0]
                     for s in np.random.SeedSequence(seed).spawn(count)],dtype=np.uint64)

def correlation(reference,image):
    reference,image=np.asarray(reference),np.asarray(image)
    valid = np.isfinite(reference)&np.isfinite(image)
    if valid.sum()<2 or np.std(reference[valid])==0 or np.std(image[valid])==0:
        return np.nan
    return float(abs(np.corrcoef(reference[valid],image[valid])[0,1]))

def summary(single,unmixed):
    """Summarize each method independently with median and linear quartiles."""
    single,unmixed=np.asarray(single),np.asarray(unmixed)
    if single.shape!=unmixed.shape or single.ndim!=2 or not len(single) or np.isinf(single).any() or np.isinf(unmixed).any():
        raise ValueError('Metrics must be paired 2D arrays with finite values or NaN')
    result = {}
    with warnings.catch_warnings():
        warnings.simplefilter('ignore',RuntimeWarning)  # all-invalid columns remain NaN
        for name,values in [('single',single),('unmixed',unmixed)]:
            result[name] = {key:np.nanquantile(values,q,axis=0,method='linear').tolist()
                            for key,q in [('median',.5),('q1',.25),('q3',.75)]}
    return result

def cache_inverse_operator(solver, measurement_count):
    """Cache a fixed linear EITLab solver's action on all measurement columns."""
    return np.column_stack([solver.solve(column) for column in np.eye(measurement_count)])


class TrialArchive:
    """Collect trial outputs, noise identities and fixed representative arrays only."""
    def __init__(self, count, names, seeds, representatives):
        self.names = np.asarray(names)
        self.seeds = seeds
        self.selected = np.array([i for i in representatives if i < count], dtype=int)
        self.single = np.empty((count, len(names)))
        self.unmixed = np.empty_like(self.single)
        self.stored = {key: [] for key in ['single_parameters', 'unmixed_parameters',
                       'single_images', 'component_images', 'unmixed_images']}
        self.noise_hashes = []

    def record(self, trial, noise, result):
        self.single[trial] = result['single_correlations']
        self.unmixed[trial] = result['unmixed_correlations']
        self.noise_hashes.append(hashlib.sha256(noise.tobytes()).hexdigest())
        if trial in self.selected:
            for key in self.stored:
                self.stored[key].append(result[key])

    def arrays(self):
        return dict(single_correlations=self.single, unmixed_correlations=self.unmixed,
                    run_seeds=self.seeds, noise_sha256=np.array(self.noise_hashes),
                    representative_trial_indices=self.selected, inclusion_names=self.names,
                    **{key: np.stack(values) for key, values in self.stored.items()})


def save_condition(directory,arrays,metadata):
    directory = Path(directory); directory.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(directory/'trials.npz',**arrays)
    metadata = dict(metadata,array_shapes={k:list(v.shape) for k,v in arrays.items()},npz_sha256=digest(directory/'trials.npz'))
    (directory/'metadata.json').write_text(json.dumps(metadata,indent=2,allow_nan=False)+'\n')
    s = summary(arrays['single_correlations'],arrays['unmixed_correlations'])
    (directory/'summary.json').write_text(json.dumps(s,indent=2)+'\n')
    return load_condition(directory)

def load_condition(directory):
    """Fail closed on inconsistent provenance, labels, axes, or stored summaries."""
    directory=Path(directory); m=json.loads((directory/'metadata.json').read_text())
    if m['schema_version']!=2 or m['algorithm']!='paired-tikhonov-v2': raise ValueError('Unknown schema/algorithm')
    if digest(directory/'trials.npz')!=m['npz_sha256']: raise ValueError('NPZ identity mismatch')
    with np.load(directory/'trials.npz',allow_pickle=False) as z: a={k:z[k] for k in z.files}
    n=m['trial_count']; p=m['parameter_count']; r=len(m['representative_trials'])
    if m['mode'] not in ('demo','full') or not isinstance(n,int) or n<1 or (m['mode']=='full' and n!=100):raise ValueError('Invalid mode/trial count')
    if m['position'] not in (.4,.6) or not np.isfinite(m['noise_level']) or m['noise_level']<0 or not np.isfinite(m['lambda_value']) or m['lambda_value']<0:raise ValueError('Invalid scientific configuration')
    if len(set(m['representative_trials']))!=r or any(not isinstance(i,int) or i<0 or i>=n for i in m['representative_trials']):raise ValueError('Invalid representative selection')
    shapes=dict(single_correlations=(n,3),unmixed_correlations=(n,3),
        run_seeds=(n,),noise_sha256=(n,),inclusion_names=(3,),representative_trial_indices=(r,),
        single_parameters=(r,p,3),unmixed_parameters=(r,p,6),single_images=(r,3,64,64),
        component_images=(r,6,64,64),unmixed_images=(r,3,64,64))
    for key,shape in shapes.items():
        if key not in a or a[key].shape!=shape or m['array_shapes'].get(key)!=list(shape): raise ValueError(f'Invalid {key}')
    if not np.array_equal(a['inclusion_names'],NAMES) or not np.array_equal(a['representative_trial_indices'],m['representative_trials']): raise ValueError('Labels/representatives mismatch')
    if not np.array_equal(a['run_seeds'],run_seeds(m['seed'],n)): raise ValueError('RNG identity mismatch')
    for key in ['single_correlations','unmixed_correlations']:
        values=a[key];finite=values[np.isfinite(values)]
        if np.isinf(values).any() or np.any((finite<0)|(finite>1+1e-12)):raise ValueError('Invalid correlation values')
    s=summary(a['single_correlations'],a['unmixed_correlations'])
    saved=json.loads((directory/'summary.json').read_text())
    for method in s:
        for key,value in s[method].items():
            if not np.allclose(value,saved[method][key],equal_nan=True,rtol=0,atol=0): raise ValueError('Summary mismatch')
    return a,m,s

def require_compatible(metadata):
    """Different position/noise may be compared; all other study settings must match."""
    keys=['algorithm','mode','trial_count','seed','lambda_value','lambda_convention','sampling','protocol','reference','basis_sha256','rng_policy','metric_policy']
    for key in keys:
        if any(m[key]!=metadata[0][key] for m in metadata[1:]): raise ValueError(f'Incompatible {key}')
