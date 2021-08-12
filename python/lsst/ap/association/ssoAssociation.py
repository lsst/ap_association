import numpy as np
import pandas as pd

from scipy.spatial import cKDTree
from .transforms import radec2icrfu


__all__ = ['associateDIA2SSO']

def associateDIA2SSO(diadf,ssodf,rOnSky=2,diara='ra',diadec='decl',ssora='AstRA(deg)',ssodec='AstDec(deg)',
                    diaId='diaSourceId',ssoId='ObjID',dId='dOnSky(arcsec)'):
    """
    Match DIA sources within rOnSky [arcsec] to Solar System Object predicted positions.
    
    Parameters:
    -----------
    diadf  ...      pandas DataFrame containing DIA source RADEC coordinates[deg]
    ssodf  ...      pandas DataFrame containing predicted SSO RADEC coordinates [deg]
    rOnSky ...      radius of kdTree search for nearest DIA neighbors to SSOs in FOV [arcsec]
    
    
    Returns:
    --------
    dia2ssodf   ... pandas DataFrame with DIA sources that have SSOs attributed to them 
                    (one SSO may have multiple dia sources in its vicinity)
    sso2dia     ... list of SSOs with attributed DIA sources
    
    
    External:
    ---------
    
    numpy, pandas, radec2icrfu, scipy: cKDTree
    
    """
    deg2rad = np.deg2rad
    rad2deg = np.rad2deg
    array = np.array
    where = np.where
    full = np.full
    norm = np.linalg.norm
    
    r=deg2rad(rOnSky/3600)
    
    # Transform DIA RADEC coordinates to unit sphere xyz for tree building.
    diaSourcesX = radec2icrfu(diadf[diara].values,diadf[diadec].values).T
    
    # Create KDTree of DIA sources
    tree = cKDTree(diaSourcesX,balanced_tree=True)
    
    treeQuery = tree.query_ball_point
    
    diaSourcesId = diadf[diaId].values
    
    sso2dia = []
    sso2diaAdd = sso2dia.append
    
    diaSsoId = []
    diaSsoIdAdd = diaSsoId.append
    
    # Query the KDtree for DIA nearest neighbors to SSOs
    for index, row in ssodf.iterrows():

        # convert SSO radec to ICRF position
        x = radec2icrfu(row[ssora],row[ssodec])

        # Which DIA Sources fall within r? 
        idx = treeQuery(x, r, p=2., eps=0)
        
        # calculate approximate on sky distance d [rad]
        d=3600.*rad2deg([norm(x-diaSourcesX[i]) for i in idx])
        
        if(len(idx)>0):
            sso2diaAdd([row[ssoId].astype(int), diaSourcesId[idx],d])
            diaSsoIdAdd(array([diaSourcesId[idx],full(len(idx),row[ssoId].astype(int)),d]).T)
        
    dia2ssodf=pd.DataFrame(np.concatenate(diaSsoId),columns=[diaId,ssoId,dId]).astype({diaId:int,ssoId:int,dId:float})

    return dia2ssodf,sso2dia
