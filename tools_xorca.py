"""Some functions with grid-aware data sets."""

import xarray as xr
import xgcm 
import matplotlib.pyplot as plt
import numpy as np

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

def plot_2D_averages(ds):
    for var in ds.data_vars:
        if len(ds[var].dims)>3:
            ax=ds[var].plot(x='t',y='depth_c',row='basins',col='data_type',col_wrap=2, sharex=False, figsize=[24,48])
            for (i,g) in enumerate(ax.axes.flat):
                g.set_title(ds.basins[int(i/2)].values)
        else:
            ax=ds[var].plot(col='basins',hue='data_type',col_wrap=3,sharey=False, sharex=False, figsize=[12,48], add_legend=False)
            for (i,g) in enumerate(ax.axes.flat):
                g.set_title(ds.basins[i].values)
                vals = g.get_yticks()
                g.set_yticklabels(['{:04.2f}'.format(x) for x in vals])
        plt.savefig(exp+'_'+var+'.png')
        
def plot_3D_averages(ds):
    for var in ds.data_vars:
        ax=ds[var].plot(row='basins',col='depth_range',hue='data_type',col_wrap=len(ds.depth_range),sharey=False, sharex=False, figsize=[12,40], add_legend=False)
        plt.savefig(exp+'_'+var+'.png')