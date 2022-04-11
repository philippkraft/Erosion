# -*- coding: utf-8 -*-
"""
A simple translation of the Eurosem model

Base equation (mass balance of the sediment):

    (D(AC)/D(t))+(D(QC)/D(x))-e(x,t)=q_s(x,t)
e=DR+DF
"""

import math

g = 9.81   # Gravity acceleration


def detachment_raindrops(rainfall, DT, PH, LD, k, P_s, z, h, PAVE):
    """
    DR(m3.s-1.m-1)=Detahment by raindrops:
    
    First term:The rainfall energy reaching the ground surface
    as direct throughfall
        
    - Rainfall: Rainfall intensity(mm.hr-1)
    - DT:Depth of direct throughfall(m) 
    - Second term:The energy of the leaf drainage
    - PH:the effective height of the plant canopy (m)
    - LD:Depth of leaf drainage received(m)
    - k:an index of the detachability of the soil (gJ-1)
    - P-s:particle density (kg m-3)
    - z:an exponent varying between 0·9 and 3·1, depending on soil texture but for
        which a value of 2·0 can be used for a wide range of conditions
    - h:the mean depth of the surface water layer (m)
    - PAVE:the fraction(between 0 and 1) of the soil surface covered by 
        non-erodible surfaces
        """    
    KE=(8.95 + (8.44 * math.log10(rainfall))) * DT + ((15.8 * (PH ** 0.5)) - 5.87) * LD
    DR=(k / P_s) * KE * (math.exp(-z * h)) * (1 - PAVE)
    return DR
 
def detachment_flow(J, w, v_s, TC, C):
    """
    returns :math:`DF(m^3 s^{-1} m^{-1})` =Soil detachment by runoff
    J:cohision of the soil(kPa)
    w:the width of flow (m)
    v_s:particle settling velocity (m.s-1)
    """
    if J < 1:
        Beta = 1
    else:
        Beta = 0.79 * math.exp(-0.85 * J)
    DF = Beta * w * v_s * (TC - C)        
    return DF


def transport_capacity_flow(d_50, u, s, U_s, y_c, h, P_s, q, n):
    """
    TODO: Need to write rill and interrill functions
    Transport capacity of the flow:
        
    - d_50:Median particle size of the soil (mm)
    - u:mean flow velocity (ms-1)
    - s:slope (%)
    - U_s:shear velocity (m.s-1)
    - y_c:modified Shields’ critical shear velocity(m.s-1),
    - P-s:particle density (kg m-3)
    - h:the mean depth of the surface water layer (m)
    - q: TODO: Discharge in m^3/s
    - n: TODO: What is n?

    """
    b=(19-(d_50/30))/10**4
    # - Pi:Modified_stream_power(g1·5 cm-2/3 s-4·5)
    Pi=((U_s*u)**1.5)/(h**(2/3))
    U_s_critical=(y_c*(P_s-1)*g*d_50)**0.5
    Pi_critical=((U_s_critical*u)**1.5)/(h**(2/3))
    
    #TC_rill=(((d_50+5)/0.32)**-0.6)*((10*u*s-0.4)**((d_50+5)/300)**0.25)
    TC_rill = 0.0
    TC_interrill=(b/(P_s * q)) * (((Pi - Pi_critical) ** (0.7/n) - 1) ** 5)
    return TC_rill + TC_interrill


# Define some values

