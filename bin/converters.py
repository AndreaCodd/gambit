from scipy.io import netcdf_file
import numpy as np
from datetime import datetime

from scipy.interpolate import griddata
def grepValuesByMask(xi, data, mask):
        """
        this grabs the values from data from entries with positive mask and interpolates it to numpy meshgrid xi
        
        """
        X=data.getX()
        x=[]
        y=[]
        z=[]
        values=[]
        for i in range(mask.getNumberOfDataPoints()):
            if mask.getTupleForDataPoint(i)[0] > 0:
                x.append(X.getTupleForDataPoint(i)[0])
                y.append(X.getTupleForDataPoint(i)[1])
                z.append(X.getTupleForDataPoint(i)[2])
                values.append(data.getTupleForDataPoint(i)[0])
        if len(xi) == 2:
            r=griddata((np.array(x), np.array(y)), np.array(values), tuple(xi), method='linear', fill_value=np.nan, rescale=False)
        else:
            r=griddata((np.array(x), np.array(y), np.array(z)), np.array(values), xi, method='linear', fill_value=np.nan, rescale=False)
        return r, xi
    

def writeNetCDF(filename, 
                data,
                error=None,
                origin=(0.,0.),
                delta=(1000.,1000.),
                units='deg',
                units_data='mgal',
                title="custom_data",
                name='data',
                longname='Data',
                summary="none",
                license="free to use",
                missing=np.nan):
    """
    create NetCDF file 
    
    :param filename: file name. include extension
    :param data: data array
    :param error: associated error. can be None, a float or an array with same size as `data`
    :param origin: tuple of origin
    :param delta: tuple of increments (can be negative)
    :param units: `deg` or `m`
    :param units_data: units of data e.g 'mgal', 'nT'
    :param title: title 
    :param name: data name
    :param longname: long data name
    :param summary: summary text
    :param license: license text
    :param missing: value for missing values
    """
    NY, NX=data.shape
    ORIGIN_X=origin[0]
    ORIGIN_Y=origin[1]
    DELTA_X=delta[0]
    DELTA_Y=delta[1]
    if units == 'm':
        XTAG='x'
        YTAG='y'
        UNITS_X="m"
        UNITS_Y="m"
    else:
        XTAG='Longitude'
        YTAG='Latitude'
        UNITS_X="degrees_east"
        UNITS_Y="degrees_north"    
    if isinstance(error, np.ndarray):
        assert error.shape == data.shape
    elif error is not None:
        print(error)
        error = np.full(data.shape, error,  dtype=data.dtype)
        error[data == missing ] = missing
        
    history=datetime.now().strftime("%d-%m-%Y")+" created using python script"

    # Create the output file and write a few metadata entries
    o=netcdf_file(filename,'w')
    o.Conventions="CF-1.0, COARDS, Unidata Dataset Discovery v1.0"
    o.Metadata_Conventions="CF-1.0, COARDS, Unidata Dataset Discovery v1.0"
    o.history=history
    o.license=license
    o.summary=summary
    o.title=title

    # Create longitude dimension and variable
    if DELTA_X > 0:
        longitude=np.linspace(ORIGIN_X, ORIGIN_X+(NX-1)*DELTA_X, NX, endpoint=True, dtype=data.dtype)
    elif DELTA_X < 0:
        longitude=np.linspace(ORIGIN_X-(NX-1)*DELTA_X, ORIGIN_X, NX, endpoint=True, dtype=data.dtype)
        
    o.createDimension(XTAG.lower(), NX)
    v=o.createVariable(XTAG.lower(), longitude.dtype, [XTAG.lower()])
    v.data[:]=longitude
    v.units=UNITS_X
    v.long_name=XTAG

    # Create latitude dimension and variable
    if DELTA_Y > 0:
        latitude=np.linspace(ORIGIN_Y, ORIGIN_Y+(NY-1)*DELTA_Y, NY, endpoint=True, dtype=data.dtype)
    elif DELTA_Y < 0:
        latitude=np.linspace(ORIGIN_Y-(NY-1)*DELTA_Y, ORIGIN_Y, NY, endpoint=True, dtype=data.dtype)
        
    o.createDimension(YTAG.lower(), NY)
    v=o.createVariable(YTAG.lower(), latitude.dtype, [YTAG.lower()])
    v.data[:]=latitude
    v.units=UNITS_Y
    v.long_name=YTAG


    # Create the main data variable
    v=o.createVariable(name, data.dtype, [YTAG.lower(), XTAG.lower()])
    v.missing_value=missing
    v.data[:]=data
    v.units=units_data
    v.long_name=longname

    # Create the error variable (can be omitted)
    if error is not None:
        v=o.createVariable(name+"_error", error.dtype, [YTAG.lower(), XTAG.lower()])
        v.missing_value=missing
        v.data[:]=error
        v.units=units_data
        v.long_name=longname+"_error"
    # Close the file
    o.close()
    return filename
    


if __name__ == "__main__":
    
    # Number of data points in longitude,latitude direction
    NX=20
    NY=10

    # Dummy value (for unset areas)
    

    # Data error (can be constant or variable over the data points)
    SIGMA = 3.

    # The actual data array, must have shape (NY, NX).
    # These are just some random numbers.
    DATA = 10*np.random.normal(size=(NY, NX), scale=SIGMA)

    # output filename
    FILENAME='test.nc'

    # Origin longitude (degrees east) and latitude (degrees north)
    ORIGIN_X=130.2
    ORIGIN_Y=-29.1

    # spacing in longitude,latitude direction (degrees)
    DELTA_X=0.05
    DELTA_Y=0.05

    # Number of data points in longitude,latitude direction
    NX=20
    NY=10

    # Data error (can be constant or variable over the data points)
    SIGMA = 3.

    # The actual data array, must have shape (NY, NX).
    # These are just some random numbers.
    DATA = 10*np.random.normal(size=(NY, NX), scale=SIGMA)

    n=writeNetCDF(filename=FILENAME, 
                data=DATA,
                units='deg',
                units_data='mgal',
                error=SIGMA,
                title='test data')
    print(f"data written to file {n}")
    

