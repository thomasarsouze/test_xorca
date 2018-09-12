"""Some functions with grid-aware data sets."""

import xarray as xr
import xgcm 

def complete_dataset(ds):
    grid = xgcm.Grid(ds, periodic=["Y", "X"])

    ds.coords['tarea'] = ds.e1t * ds.e2t
    ds.coords['uarea'] = ds.e1u * ds.e2u
    ds.coords['varea'] = ds.e1v * ds.e2v
    ds.coords['farea'] = ds.e1f * ds.e2f
    if 'thkcello' in ds.variables.keys():
        ds.coords['e3t'] = ds.thkcello
        ds.coords['e3u'] = grid.interp(ds.thkcello,"X",boundary="fill")
        ds.coords['e3v'] = grid.interp(ds.thkcello,"Y",boundary="fill")
        ds.coords['e3w'] = grid.interp(ds.thkcello,"Z",boundary="fill")
    ds.coords['tvol'] = ds.tarea * ds.e3t
    ds.coords['uvol'] = ds.uarea * ds.e3u
    ds.coords['vvol'] = ds.varea * ds.e3v
    ds.coords['wvol'] = ds.tarea * ds.e3w

    return ds
