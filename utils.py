import constants as cnst
import math

def getSurfaceVelocity(I, A=None)
    I *= 1E-09 # in A from nA
    A = cnst.A_GCD if A is None else A
    A *= 1E-6 # in m^2
    ni = cnst.ni * 1E+06 # in m^-3
    J = I / (cnst.q*ni*A) # in m/s
    return J*100 # in cm/s

def getOxideChargeDensity(V_fb, structure, C_ox, A=None, approx_openC=50): 
    if structure not in ['MOShalf', 'MOS2000']:
        raise ValueError(f"Invalid structure: {structure}. Expected 'MOShalf' or 'MOS2000'.")
    if A is None:
        A = cnst.A_GCD if structure == 'MOShalf' else cnst.A_MOS2000
    A *= 1E-6 
    C_ox -= approx_openC # in pF
    C_ox *= 1E-12 # in F from pF
    phi_s = cnst.Chi + cnst.Eg/2. + cnst.KbTn20 * math.log(cnst.NA/cnst.ni)
    phi_ms = cnst.phi_m - phi_s
    N_ox = C_ox / (A*cnst.q) * (phi_ms + V_fb) # 1/m^2
    return N_ox * 1E-04 # in 1/cm^2


