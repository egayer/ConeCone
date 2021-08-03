#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  3 18:24:13 2021

@author: egayer
"""

import numpy as np
from scipy.optimize import fmin
from scipy.stats import spearmanr


#%% Compute distance between 2 points

def radius(x,y,z,xa,ya) :
    """ 
    Compute distances between the control points (x,y) and a potential center (xa, ya) 
    in UTM coord, and sort the results depending on the elevation 
    of the control points (z)
    
    Parameters
    ----------
    x,y,z : UTM coordinates and Elevation of the control point
    xa, ya : UTM coordinates of the potential center ("a" for apex)

    Returns
    -------
    r: Distance between the control point and the potential center 
    ("r" for Radius since we are dealing with conical shapes)
    z: Elevation of the control point
    Both r and z are sorted depending on the elevation.
        
    Example
    -------
    radius_sorted, Elev_sorted = radius(control_xUTM, control_yUTM,control_z,apex_coordX,apex_coordY)
    
    """
    r = np.sqrt( (x-xa)**2 + (y-ya)**2 )
    index = np.argsort(r)
    return r[index], np.array(z)[index]


#%% Test Spearman correlation

def errorSpearman(r,z_sorted) :
    """
    Calculate Spearman's coefficient rho
    between the distances (center<->control points)
    and the elevation of the control points
    
    Parameters
    ----------
    r : Distance or radius between the control points and the potential center (sorted)
    z_sorted : Elevation of the control points (sorted)

    Returns
    -------
    rho : Spearman's coefficient rho
    
    """
    rho, pval = spearmanr(r,z_sorted)
    return rho



#%% Test Linear Regression

def errorPolynome(r,z_sorted,deg=1) :
    """
    Calculate the errors of a regression line (when deg=1) fitting 
    the distances (center<->control points)
    and the elevation of the control points
    
    Parameters
    ----------
    r : Distance or radius between the control points and the potential center (sorted)
    z_sorted : Elevation of the control points (sorted)

    Returns
    -------
    residuals
    """
    p, residuals, rank, singular_values, rcond = np.polyfit(r,z_sorted,deg,full = True)
    #return sum(residuals**2)
    return residuals


#%%

def conecone_CenterProfile(control_xUTM_m,control_yUTM_m, control_z_m, x_center_guess, y_center_guess):
    """
    Gives the best center and the corresponding distance between the best
    center and every control points, for both 
    -1- LR (Linear Regression) method which uses a 1 degree polynomial to fit the control points
    -2- Sp (Spearman Correlation) method which uses the spearman rank correlation to estimate the 
    best ordering of the control points
    
    Parameters
    ----------
    control_xUTM_m, control_yUTM_m, control_z_m : UTM coordinates and Elevation of the control points
    x_center_guess, y_center_guess : UTM coordinates of the potential center

    Returns
    -------
    Results are given in a Dictionnary with the (x,y) of the best centers of the structure which are
    the locations that give the both the best Linear Regression coeficient and the bes Spearman coefficient
    The dictionary contains also the distance to the best centers and the elevations of all the
    control points.
    
    """
    
    def to_be_minimizedP1( apex_coord ) :
        rP1, zrP1 = radius(control_xUTM_m, control_yUTM_m,control_z_m,apex_coord[0],apex_coord[1])
        return errorPolynome(rP1,zrP1,1)

    def to_be_minimizedSp( apex_coord ) :
        rSp, zrSp = radius(control_xUTM_m, control_yUTM_m,control_z_m,apex_coord[0],apex_coord[1])
        return errorSpearman(rSp,zrSp)

    BestCenter_P1 = fmin(func=to_be_minimizedP1, x0=[x_center_guess,y_center_guess],disp=False)
    BestCenter_Sp = fmin(func=to_be_minimizedSp, x0=[x_center_guess,y_center_guess],disp=False)
    
    r_P1, zr_P1 = radius(control_xUTM_m, control_yUTM_m,control_z_m,BestCenter_P1[0],BestCenter_P1[1])
    r_Sp, zr_Sp = radius(control_xUTM_m, control_yUTM_m,control_z_m,BestCenter_Sp[0],BestCenter_Sp[1])
    
    resDic = dict()  
    resDic['BestCenter_LR_xy'] = BestCenter_P1
    resDic['RadiusModel_controlPoint_LR'] = r_P1
    resDic['Elevation_controlPoint_LR'] = zr_P1
    
    resDic['BestCenter_Sp_xy'] = BestCenter_Sp
    resDic['RadiusModel_controlPoint_Sp'] = r_Sp
    resDic['Elevation_controlPoint_Sp'] = zr_Sp
    
    return resDic