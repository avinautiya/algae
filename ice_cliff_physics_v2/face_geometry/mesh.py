"""Lossless insertion of supplied class interfaces in an extruded face profile.

No coarsening, heat transfer, ice removal or material transport is performed.
The maximum-band contract is checked before constructing a replacement profile.
"""
from bisect import bisect_left, bisect_right
from dataclasses import replace
try:
    from .face_profile import Band, Profile, finite
except ImportError:
    from face_profile import Band, Profile, finite


def split_at_elevations(profile, elevations, max_bands=4096):
    """Return the identical geometric body on a possibly refined height mesh.

    Sorted unique finite elevations strictly inside an existing band become new
    interfaces. Both child fronts copy their parent's coordinate exactly. Values
    outside the face and existing interfaces do not add bands. No epsilon-based
    merging/coarsening is applied, so callers can audit very close interfaces.
    The input Profile and the caller's elevations are never modified.
    """
    if not isinstance(profile, Profile):
        raise ValueError('validated Profile required')
    if isinstance(max_bands, bool) or not isinstance(max_bands, int) or max_bands <= 0:
        raise ValueError('max_bands must be a positive integer')
    if isinstance(elevations, (str, bytes, dict)):
        raise ValueError('finite elevation iterable required')
    try:
        points=sorted(set(finite(z,'split elevation') for z in elevations))
    except TypeError as exc:
        raise ValueError('finite elevation iterable required') from exc
    slices=[];count=len(profile.bands)
    for b in profile.bands:
        lo=bisect_right(points,b.z_low_m)
        hi=bisect_left(points,b.z_high_m)
        slices.append((lo,hi));count+=hi-lo
    if count>max_bands:
        raise ValueError(f'height mesh would require {count} bands, exceeding cap {max_bands}; no coarsening or mutation performed')
    if count==len(profile.bands):
        return profile
    bands=[]
    for b,(lo,hi) in zip(profile.bands,slices):
        if lo==hi:
            bands.append(b)
            continue
        edges=[b.z_low_m]+points[lo:hi]+[b.z_high_m]
        bands.extend(Band(a,z,b.front_x_m) for a,z in zip(edges,edges[1:]))
    return replace(profile,bands=tuple(bands))
