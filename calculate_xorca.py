"""Calculations with grid-aware data sets."""

import xgcm
from __param__ import *
import xarray as xr

def calculate_moc(grid,ds, region=""):
    """Calculate the MOC.
    Parameters
    ----------
    grid : grid associated with ds
    ds : xarray dataset
        A grid-aware dataset as produced by `xorca.lib.preprocess_orca`.
    region : str
        A region string.  Examples: `"atl"`, `"pac"`, `"ind"`.
        Defaults to `""`.
    Returns
    -------
    moc : xarray data array
        A grid-aware data array with the moc for the specified region.  The
        data array will have a coordinate called `"lat_moc{region}"` which is
        the weighted horizontal and vertical avarage of the latitude of all
        latitudes for the given point on the y-axis.
    """
    vmaskname = "vmask" + region
    mocname = "moc" + region
    latname = "lat_moc" + region

    weights = ds[vmaskname] * ds.e3v * ds.e1v

    Ve3 = weights * ds.vo

    # calculate indefinite vertical integral of V from bottom to top, then
    # integrate zonally, convert to [Sv], and rename to region
    moc = grid.cumsum(Ve3, "Z", to="left", boundary="fill") - Ve3.sum("z_c")
    moc = moc.sum("x_c")
    moc /= 1.0e6
    moc = moc.rename(mocname)

    # calculate the weighted zonal and vertical mean of latitude
    lat_moc = ((weights * ds.llat_rc).sum(dim=["z_c", "x_c"]) /
               (weights).sum(dim=["z_c", "x_c"]))
    moc.coords[latname] = (["y_r", ], lat_moc.data)

    # also copy the relevant depth-coordinates
    moc.coords["depth_l"] = ds.coords["depth_l"]
    
    moc = moc.compute()

    return moc


def calculate_psi(grid,ds):
    """Calculate the barotropic stream function.
    Parameters
    ----------
    grid : grid associated with ds
    ds : xarray dataset
        A grid-aware dataset as produced by `xorca.lib.preprocess_orca`.
    Returns
    -------
    psi : xarray data array
        A grid-aware data array with the barotropic stream function in `[Sv]`.
    """
    U_bt = (ds.uo * ds.e3u).sum("z_c")

    psi = grid.cumsum(- U_bt * ds.e2u, "Y") / 1.0e6
    psi -= psi.isel(y_r=-1, x_r=-1)  # normalize upper right corner
    psi = psi.rename("psi")

    psi = psi.compute()

    return psi


def calculate_speed(grid,ds):
    """Calculate speed on the central (T) grid.
    First, interpolate U and V to the central grid, then square, add, and take
    root.
    Parameters
    ----------
    grid : grid associated with ds
    ds : xarray dataset
        A grid-aware dataset as produced by `xorca.lib.preprocess_orca`.
    Returns
    -------
    speed : xarray data array
        A grid-aware data array with the speed in `[m/s]`.
    """
    U_cc = grid.interp(ds.uo, "X", to="center")
    V_cc = grid.interp(ds.vo, "Y", to="center")

    speed = (U_cc**2 + V_cc**2)**0.5

    speed = speed.compute()

    return speed

def calculate_ke(grid,ds,full=False,depths=[0]):
    """Calculate KE, MKE, EKE on the central (T) grid.
    First, interpolate U and V to the central grid, then square, add, and multiply by 0.5.
    Parameters
    ----------
    grid : grid associated with ds
    ds : xarray dataset
        A grid-aware dataset as produced by `xorca.lib.preprocess_orca`.
    full : by default, only MKE and EKE in the output. If set to true, also put KE at each time step
    levels : list of depths over which to calculate the KE budget. By default only the surface.
    Returns
    -------
    speed : xarray dataset
        A grid-aware data array with the speed in `[m/s]`.
        Includes KE (if requested)
    """
    arrays=[]
    for depth in depths:
        U_d = ds.uo.sel(depth_c=depth, method='nearest')
        V_d = ds.vo.sel(depth_c=depth, method='nearest')
        U_cc = grid.interp(U_d, "X", to="center")
        V_cc = grid.interp(V_d, "Y", to="center")
        U_cc_m = U_cc.mean(dim='t')
        V_cc_m = V_cc.mean(dim='t')
        if 'u2o' in ds.variables:
            U2_d = ds.u2o.sel(depth_c=depth, method='nearest')
            V2_d = ds.v2o.sel(depth_c=depth, method='nearest')
            U2_cc = grid.interp(U_d, "X", to="center")
            V2_cc = grid.interp(V_d, "Y", to="center")
            TKE = (0.5*(U2_cc + V2_cc)).mean(dim='t')
            MKE = 0.5*(U_cc_m**2 + V_cc_m**2)
            EKE = TKE - MKE
        else:
            MKE = 0.5*(U_cc_m**2 + V_cc_m**2)
            EKE = (0.5*((U_cc-U_cc_m)**2 + (V_cc-V_cc_m)**2)).mean(dim='t')
        MKE.name='MKE_'+str(depth)+'m'
        EKE.name='EKE_'+str(depth)+'m'
        MKE = MKE.compute()
        EKE = EKE.compute()
        arrays.append(MKE)
        arrays.append(EKE)
        if full:
            KE  = 0.5*(U_ccd**2 + V_ccd**2)
            KE.name='KE_'+str(depth)+'m'
            KE = KE.compute()
            arrays.append(KE)

    return xr.merge(arrays)

def calculate_enso(da):
    """Calculate ENSO (over Nino3.4 box)
    First select the zone, then calculate the mean SST, then the SST anomaly using a 5 months filter
    Parameters
    ----------
    ds : xarray dataarray of SST
        A grid-aware dataset as produced by `xorca.lib.preprocess_orca`.
    Returns
    -------
    enso : xarray data array
        A grid-aware data array with the time series of the indice
    """
    lon1_nino = - 170. ; lon2_nino = - 120. ; lat1_nino = -5. ; lat2_nino = 5.
    condition = ((da.llon_cc>lon1_nino) & (da.llon_cc<lon2_nino) & (da.llat_cc>lat1_nino) & (da.llat_cc<lat2_nino)).squeeze()
    da = da.where(condition,drop=True)
    dims = da.dims[-2:]
    sst  = (v * v['tarea'] ).sum(dims) / (v['area']).sum(dims)
    sst_anom = sst - sst.mean(dim='t')
    enso = sst_anom.rolling(t=5, center=True).mean()

    enso = enso.compute()
    
    return enso

def calculate_wind_stress_curl(grid,ds):
    """Calculate wind-stress curl on vorticity F grid.
    Calculate the curl (dTy/dx - dTx/dy)
    Parameters
    ----------
    grid : grid associated with ds
    ds : xarray dataset
        A grid-aware dataset as produced by `xorca.lib.preprocess_orca`.
    Returns
    -------
    curl : xarray dataset
        A grid-aware data array with the wind stress curl in `[N/m3]`.
    """

    curl = (grid.diff(ds.tauvo * ds.e1v, 'X') + grid.diff(ds.tauuo * ds.e2u, 'Y'))/ds.farea
    curl = curl.compute()

    return curl

def calculate_vorticity(grid,ds):
    """Calculate vorticity on vorticity F grid.
    Calculate the strain (dV/dx - dU/dy)
    Parameters
    ----------
    grid : grid associated with ds
    ds : xarray dataset
        A grid-aware dataset as produced by `xorca.lib.preprocess_orca`.
    levels : list of depths over which to calculate the KE budget. By default only the surface.
    Returns
    -------
    curl : xarray dataset
        A grid-aware data array with the wind stress curl in `[N/m3]`.
    """

    vorticity = (grid.diff(ds.vo * ds.e1v, 'X') + grid.diff(ds.uo * ds.e2u, 'Y'))/ds.farea
    vorticity = vorticity.compute()

    return vorticity

def calculate_strain(grid,ds):
    """Calculate strain on central T grid.
    Calculate the vorticity (dU/dx - dV/dy)
    Parameters
    ----------
    grid : grid associated with ds
    ds : xarray dataset
        A grid-aware dataset as produced by `xorca.lib.preprocess_orca`.
    levels : list of depths over which to calculate the KE budget. By default only the surface.
    Returns
    -------
    curl : xarray dataset
        A grid-aware data array with the wind stress curl in `[N/m3]`.
    """

    strain = (grid.diff(ds.uo * ds.e2u, 'X') + grid.diff(ds.vo * ds.e1v, 'Y'))/ds.tarea
    strain = strain.compute()

    return strain


def _read_sections(sections):
    f = open(sections, "r") 
    name=[];indices=[];ref_TS=[];
    while True:
        name.append(f.readline())
        if 'EOF' in name[-1]: break
        indices.append(f.readline().split())
        last_pos = f.tell()
        if "ref_temp_sali" in f.readline():
            f.seek(last_pos)
            ref_TS.append(f.readline().split(":")[-1].split())
        else:
            f.seek(last_pos)
            ref_TS.append([0.,34.8])
            name = f.readline()
    
    f.close()
    
    return name,indices,ref_TS

def _get_direction(indice):
    i1,i2,j1,j2=indice
    idirx = np.sign(i2-i1)
    idirj = np.sign(j2-j1)
    norm_u =  idiry
    norm_v = -idirx

    return norm_u, norm_v

def _get_points(section,dims):
    """
    Reads in a file the list of points to be considered for the current section. The file to be read is obtained from cdftool cdftransportiz. Needs also dims : the disctionnary of dimensions from the dataset.
    Outputs mask_u and mask_v
    """
    fname="/home/Earth/tarsouze/tmp/broken_line_"+section+".dat"
    points=[]
    maskMFO_u = xr.DataArray(np.zeros((dims['y_c'],dims['x_r']), dtype=np.int16),dims=('y_c','x_r'),name= 'maskMFO_u')
    maskMFO_v = xr.DataArray(np.zeros((dims['y_r'],dims['x_c']), dtype=np.int16),dims=('y_r','x_c'),name= 'maskMFO_v')
    with open(fname) as f:
        next(f)
        for line in f:
            points.append(line.split())

    for l in range(len(points)-1):
        if (points[l+1][0]==points[l][0]):
           if (int(points[l+1][1])==int(points[l][1])+1):
               maskMFO_u[int(points[l][1]),int(points[l][0])-1]=1
           else:
               maskMFO_u[int(points[l][1])-1,int(points[l][0])-1]=1
        if (points[l+1][1]==points[l][1]):
           if (int(points[l+1][0])==int(points[l][0])+1):
               maskMFO_v[int(points[l][1])-1,int(points[l][0])]=1
           else:
               maskMFO_v[int(points[l][1])-1,int(points[l][0])-1]=1

    return masksMFO_u, maskMFO_v



def calculate_transport_sections(ds,sections,depths=[0]):
    """Calculate transport through sections
    First read the definition of the sections, then calculate volume transport, heat and salt transport, finaly extract transport on section points.
    Parameters
    ----------
    ds : xarray dataset
        A grid-aware dataset as produced by `xorca.lib.preprocess_orca`.
    sections : name of the files describing the sections. Same file as for cdftools.
    depths : list of depths defining depth layers for calculation
    Returns
    -------
    transport : xarray dataset
        A grid-aware dataset with the time series of volume, heat, salt transport for each section
    """
    depths.append(10000)
    depth = [key for key in ds.coords.keys() if key.startswith('depth')][0]

    # Compute the transport and sum
    names,indices,refs_TS = _read_sections(sections)
    trpu = ds.uo * ds.e3u * ds.e2u
    trpv = ds.vo * ds.e3v * ds.e1v
    trput = trpu * grid.interp(ds_xorca.thetao,'X') * rau0 * rcp
    trpvt = trpv * grid.interp(ds_xorca.thetao,'Y') * rau0 * rcp
    trpus = trpu * grid.interp(ds_xorca.so,'X')
    trpvs = trpv * grid.interp(ds_xorca.so,'Y')
    array_all=[]
    # loop over sections
    for name,indice,ref_TS in zip(names,indices,refs_TS):
        mask_u, mask_v = _get_points(name,ds.dims)
        norm_u, norm_v = _get_direction(indice)
        # then loop over depths
        if len(depths)>2:
            for depth1,depth2 in zip(depths,depths[1:]):
                condition = ((-ds[depth]>depth1) & (-ds[depth]<depth2))
                arrayv = trpu.where(condition & mask_u).sum(trpu.dims[-3:])*norm_u + trpv.where(condition & mask_v).sum(trpv.dims[-3:])*norm_v 
                arrayv = arrayv.compute()
                arrayt = trput.where(condition & mask_u).sum(trput.dims[-3:])*norm_u + trpvt.where(condition & mask_v).sum(trpvt.dims[-3:])*norm_v
                arrayt -= rau0 * rcp * ref_TS[0] * arrayv 
                arrayt = arrayt.compute()
                arrays = trpus.where(condition & mask_u).sum(trpus.dims[-3:])*norm_u + trpvs.where(condition & mask_v).sum(trpvs.dims[-3:])*norm_v
                arrayf = arrayv - arrays/ref_TS[-1] # only valid if ref_sali=0.
                arrays -= ref_TS[-1] * arrayv
                arrays = arrays.compute()
                if depth2==10000:
                    arrayv.name='trp_volume_'+name+'_'+str(depth1)+'-bottom'
                    arrayt.name='trp_heat_'+name+'_'+str(depth1)+'-bottom'
                    arrays.name='trp_salt_'+name+'_'+str(depth1)+'-bottom'
                    arrayf.name='trp_freshwater_'+name+'_'+str(depth1)+'-bottom'
                else:
                    arrayv.name='trp_volume_'+name+'_'+str(depth1)+'-'+str(depth2)
                    arrayt.name='trp_heat_'+name+'_'+str(depth1)+'-'+str(depth2)
                    arrays.name='trp_salt_'+name+'_'+str(depth1)+'-'+str(depth2)
                    arrayf.name='trp_freshwater_'+name+'_'+str(depth1)+'-'+str(depth2)
                arrayv.attrs['units']='Sv'
                arrayt.attrs['units']='PW'
                arrayt.attrs['Tref']=ref_TS[0]
                arrays.attrs['units']='kt/s'
                arrays.attrs['Sref']=ref_TS[-1]
                arrayf.attrs['units']='Sv'
                arrayf.attrs['Sref']=ref_TS[-1]
                array_all.append((arrayv*1.e-6,arrayt*1.e-15,arrays*1e-6,arrayf*1e-6))

        #add the computation for the whole column too
        arrayv = trpu.where(mask_u).sum(trpu.dims[-3:])*norm_u + trpv.where(mask_v).sum(trpv.dims[-3:])*norm_v
        arrayt = trput.where(mask_u).sum(trput.dims[-3:])*norm_u + trpvt.where(mask_v).sum(trpvt.dims[-3:])*norm_v
        arrays = trpus.where(mask_u).sum(trpus.dims[-3:])*norm_u + trpvs.where(mask_v).sum(trpvs.dims[-3:])*norm_v
        arrayv.name = 'trp_volume_'+name+'_0-bottom'
        arrayt.name = 'trp_heat_'+name+'_0-bottom'
        arrays.name = 'trp_salt_'+name+'_0-bottom'
        arrayv.attrs['units']='Sv'
        arrayt.attrs['units']='PW'
        arrayt.attrs['Tref']=ref_TS[0]
        arrays.attrs['units']='kt/s'
        arrays.attrs['Tref']=ref_TS[-1]
        arrayf.attrs['units']='Sv'
        arrayf.attrs['Sref']=ref_TS[-1]
        array_all.append((arrayv*1.e-6,arrayt*1.e-15,arrays*1e-6,arrayf*1e-6))

    return xr.merge(array_all)


def average_2D(ds,var):
    """Calculate the area weighted average of 'var' variable in ds
    Parameters
    ----------
    ds : xarray dataset
        A grid-aware dataset as produced by `xorca.lib.preprocess_orca`.
    var : name of the variable to average. Has to be included in ds.
    Returns
    -------
    ave : xarray data array
        A grid-aware data array with the averaged variable. 
        Time_series format if var is 2D+time
        Time_series + depth if var is 3D+time
    """
    
    if var in ds.variables.keys():
        v = ds[var]
        dims = orca_variables[var]['dims'][-2:]
        area = [key for key in v.coords.keys() if key.endswith('area')][0]
        ave = (v * v[area] ).sum(dims) / (v[area]).sum(dims)
        ave = ave.compute()
    else:
        print('Variable '+var+' is not in the dataset. Impossible to do average_2D. Please select another variable')

    return ave

def average_3D(ds,var,depths=[0]):
    """Calculate the volume weighted average of 'var' variable in ds
    Parameters
    ----------
    ds : xarray dataset
        A grid-aware dataset as produced by `xorca.lib.preprocess_orca`.
    var : name of the variable to average. Has to be included in ds and be 3D (+ time)
    depths : list of depths defining depth layers for averaging
    Returns
    -------
    ave : xarray data array
        A grid-aware data array with the averaged variable. 
        Time_series format 
    """

    if var in ds.variables.keys():
        v = ds[var]
        depths.append(10000) # do the calculation from last depth to the bottom
        dims  = update_orca_variables[var]['dims'][-3:]
        depth = [key for key in v.coords.keys() if key.startswith('depth')][0]
        vol   = [key for key in v.coords.keys() if key.endswith('vol')][0]
        arrays=[]
        if len(depths)>2:
            for depth1,depth2 in zip(depths,depths[1:]):
                condition = ((-v[depth]>depth1) & (-v[depth]<depth2))
                array = (v * v[vol] ).where(condition,drop=True).sum(dims) / (v[vol]).where(condition,drop=True).sum(dims)
                array = array.compute()
                if depth2==10000:
                    array.name=var+'_'+str(depth1)+'-bottom'
                else:
                    array.name=var+'_'+str(depth1)+'-'+str(depth2)
                arrays.append(array)
        #add the computation for the whole column too
        array = (v * v[vol] ).sum(dims) / (v[vol]).sum(dims)
        array.name = var+'_0-bottom'
        arrays.append(array)
    else:
        print('Variable '+var+' is not in the dataset. Impossible to do average_3D. Please select another variable')

    return xr.merge(arrays) 

