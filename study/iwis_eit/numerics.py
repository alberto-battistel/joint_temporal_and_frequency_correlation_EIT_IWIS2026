"""Image sampling details; scientific reconstruction steps live in scripts 01 and 02."""
import numpy as np


def endpoint_image(field, mesh):
    """64x64 endpoint grid; preserve native lookup and finite mask."""
    image = np.full((64,64),np.nan)
    for j,y in enumerate(np.linspace(-1,1,64)):
        for i,x in enumerate(np.linspace(-1,1,64)):
            try:
                image[j,i] = field(mesh(x,y))
            except Exception:  # Native outside-mesh lookup raises a version-specific exception.
                pass
    return image



def cached_parameter_images(parameters, pixel_dofs):
    """Apply the saved native pixel lookup, with no new FEM evaluations."""
    images = np.full((parameters.shape[1], 64, 64), np.nan)
    valid = pixel_dofs >= 0
    images[:, valid] = parameters[pixel_dofs[valid], :].T
    return images
