# =============================================================================
# Copyright and citation notice
# =============================================================================
#
# This code implements the PLS-based cutoff latitude and proton cutoff energy
# model developed in the following manuscript:
#
#     A New Model of Cutoff Latitude for Solar Energetic Protons Based on the
#     Partial Least Squares Method
#
# Authors:
#     Zhenghao She^1, Hui Zhu^1,*, Qian Song^2, Yao Chen^1,*,
#     Yingying Zhao^1, Huicong Chen^1, and Zhijie Qin^1
#
# Affiliations:
#     ^1 Institute of Frontier and Interdisciplinary Science,
#        Shandong University, Shandong, China
#
#     ^2 Key Laboratory of Space Weather, National Satellite Meteorological
#        Center (National Center for Space Weather), China Meteorological
#        Administration, Beijing, China
#
# Corresponding authors:
#     * Hui Zhu, huizhu@email.sdu.edu.cn
#     * Yao Chen, yaochen@sdu.edu.cn
#
# The cutoff latitude and proton cutoff energy models implemented in this file
# are derived from the above manuscript. If you use this model in research work
# intended for publication, we recommend citing the manuscript. The final
# bibliographic information should be updated after the manuscript is published.
#
# For questions, bug reports, or problems encountered when using the model,
# please contact:
#
#     Hui Zhu, huizhu@email.sdu.edu.cn
#     Zhenghao She, shezhenghao@qq.com
#
# Copyright (c) 2026 Zhenghao She, Hui Zhu, Qian Song, Yao Chen,
# Yingying Zhao, Huicong Chen, and Zhijie Qin. All rights reserved.
#
# =============================================================================

"""
PLS-based cutoff latitude and proton cutoff energy model.

This module provides a PLS-based empirical model for calculating cutoff
latitude and proton cutoff energy using solar wind parameters, geomagnetic
activity indices, proton energy, geomagnetic latitude, and magnetic local
time (MLT).

The four main user-facing functions are:

```
calculate_cutoff_latitude_divisional_mlt
calculate_cutoff_latitude_continuous_mlt
calculate_proton_cutoff_energy_divisional_mlt
calculate_proton_cutoff_energy_continuous_mlt
```

1. Cutoff latitude calculation

---

The following two functions calculate cutoff latitude:

```
calculate_cutoff_latitude_divisional_mlt(
    Kp, Dst, Pdyn, SYM_H, Bt, Vt, AE, proton_energy, MLT
)

calculate_cutoff_latitude_continuous_mlt(
    Kp, Dst, Pdyn, SYM_H, Bt, Vt, AE, proton_energy, MLT
)
```

The inputs must be one-dimensional NumPy arrays of floating-point values
with the same length. Elements with the same index in all input arrays are
assumed to correspond to the same time point.

Input variables:

```
Kp
    Kp index, in decimal representation (i.e., Kp × 10).

Dst
    Dst index, in nT.

Pdyn
    Solar wind dynamic pressure, in nPa.

SYM_H
    SYM-H index, in nT.

Bt
    Total interplanetary magnetic field strength, in nT.

Vt
    Solar wind speed, in km/s.

AE
    AE index, in nT.

proton_energy
    Proton energy, in MeV. The value must be positive because the model
    uses the natural logarithm of proton energy, i.e. np.log(proton_energy).

MLT
    Magnetic local time, in hours. The valid range is 0 <= MLT <= 24.
```

The output of calculate_cutoff_latitude_divisional_mlt is:

```
cutoff_latitude_divisional_mlt
```

which is a one-dimensional NumPy array of floating-point values. It gives the
cutoff latitude calculated using the divisional MLT model. The unit is degree.
The output array has the same length as the input arrays, and each element
corresponds to the same time point as the elements with the same index in the
input arrays.

The output of calculate_cutoff_latitude_continuous_mlt is:

```
cutoff_latitude_continuous_mlt
```

which is a one-dimensional NumPy array of floating-point values. It gives the
cutoff latitude calculated using the continuous MLT model obtained by linear
interpolation in MLT. The unit is degree. The output array has the same length
as the input arrays, and each element corresponds to the same time point as
the elements with the same index in the input arrays.

2. Proton cutoff energy calculation

---

The following two functions calculate proton cutoff energy by inverting the
corresponding cutoff latitude formulas:

```
calculate_proton_cutoff_energy_divisional_mlt(
    Kp, Dst, Pdyn, SYM_H, Bt, Vt, AE, Latitude, MLT
)

calculate_proton_cutoff_energy_continuous_mlt(
    Kp, Dst, Pdyn, SYM_H, Bt, Vt, AE, Latitude, MLT
)
```

The inputs must be one-dimensional NumPy arrays of floating-point values
with the same length. Elements with the same index in all input arrays are
assumed to correspond to the same time point.

Input variables:

```
Kp
    Kp index, in decimal representation (i.e., Kp × 10).

Dst
    Dst index, in nT.

Pdyn
    Solar wind dynamic pressure, in nPa.

SYM_H
    SYM-H index, in nT.

Bt
    Total interplanetary magnetic field strength, in nT.

Vt
    Solar wind speed, in km/s.

AE
    AE index, in nT.

Latitude
    Geomagnetic latitude, in degree.

MLT
    Magnetic local time, in hours. The valid range is 0 <= MLT <= 24.
```

The output of calculate_proton_cutoff_energy_divisional_mlt is:

```
proton_cutoff_energy_divisional_mlt
```

which is a one-dimensional NumPy array of floating-point values. It gives the
proton cutoff energy calculated by inverting the divisional MLT cutoff
latitude model. The unit is MeV. The output array has the same length as the
input arrays, and each element corresponds to the same time point as the
elements with the same index in the input arrays.

The output of calculate_proton_cutoff_energy_continuous_mlt is:

```
proton_cutoff_energy_continuous_mlt
```

which is a one-dimensional NumPy array of floating-point values. It gives the
proton cutoff energy calculated by inverting the continuous MLT cutoff
latitude model. The unit is MeV. The output array has the same length as the
input arrays, and each element corresponds to the same time point as the
elements with the same index in the input arrays.

3. MLT treatment

---

The divisional MLT model uses four MLT sectors:

```
dayside:   9  <= MLT < 15
duskside: 15 <= MLT < 21
nightside: 21 <= MLT <= 24 and 0 <= MLT < 3
dawnside: 3  <= MLT < 9
```

The continuous MLT model treats the sectoral cutoff latitudes as values at
four representative MLT locations:

```
dawnside   -> 06 MLT
dayside    -> 12 MLT
duskside   -> 18 MLT
nightside  -> 00/24 MLT
```

The cutoff latitude is then obtained by linear interpolation in MLT between
these representative locations. The continuous proton cutoff energy is
obtained by inverting this continuous cutoff latitude formulation.

4. Notes

---

All calculations are performed element by element. Therefore, all input arrays
must have the same length and must be aligned in time.

Invalid or out-of-range values may produce np.nan in the output. In particular,
MLT values outside the range 0 <= MLT <= 24 are not assigned to any MLT sector.
"""


import numpy as np

################################################################

PLS1_DAYSIDE_COEFFS = np.array([-0.010029861920951807, 0.0069841242900215645, -0.19788310507019724,
                               0.0070130597086824403, -0.017521284700369961, -0.0016525262927944821,
                               -0.00054572448329867141, -0.46975213640591373], dtype=np.float64)

PLS1_DAYSIDE_INTERCEPT = np.float64(3.0971194504776061)

PLS1_DUSKSIDE_COEFFS = np.array([-0.021224437849494186, 0.0081775223314167322, -0.65571323761240385,
                                0.0079596027589553184, -0.049954345370450827, -0.0028772586992722661,
                                -0.00091238332764227552, -0.17160036123002120], dtype=np.float64)

PLS1_DUSKSIDE_INTERCEPT = np.float64(4.3489522140435497)

PLS1_NIGHTSIDE_COEFFS = np.array([-0.019777890049169009, 0.0064065471432762661, -0.74318370330909123,
                                 0.0062113399135833832, -0.048206363667649514, -0.0028923243466119604,
                                 -0.00086175141593595702, -0.090802964785902759], dtype=np.float64)

PLS1_NIGHTSIDE_INTERCEPT = np.float64(4.3408013028081198)

PLS1_DAWNSIDE_COEFFS = np.array([-0.015542987823886420, 0.0084059113554201988, -0.49791473761336169,
                                0.0076519532324660052, -0.038555100940759000, -0.0025471493230617070,
                                -0.00088079787969121972, -0.36262852322948852], dtype=np.float64)

PLS1_DAWNSIDE_INTERCEPT = np.float64(4.2338669401962843)

##################################################################

CUTOFF_LATITUDE_DAYSIDE_LEFT_COEFFS = np.array([-0.0027323584797300324, 0.011900738674032976, 0.47951755900485893,
                                               0.012101964591507530, 0.017933035625206737, -0.00057110270156969750,
                                               -0.00050800386764092030, -1.5381707084554828], dtype=np.float64)

CUTOFF_LATITUDE_DAYSIDE_LEFT_INTERCEPT = np.float64(69.273256888492980)

CUTOFF_LATITUDE_DAYSIDE_RIGHT_COEFFS = np.array([-0.017780026023950641, 0.022378926831760584, 0.18263618604256493,
                                                0.022623564170418525, -0.0083539132078621486, -0.0030503657719924008,
                                                -0.0013267470029044110, -2.2429335489146824], dtype=np.float64)

CUTOFF_LATITUDE_DAYSIDE_RIGHT_INTERCEPT = np.float64(74.129578110918203)

CUTOFF_LATITUDE_DAYSIDE_BREAKPOINT = np.float64(-0.13980953945215768)

####################################################################

CUTOFF_LATITUDE_DUSKSIDE_LEFT_COEFFS = np.array([-0.016598095069219932, 0.0095013284270082168, -0.25225924451919618,
                                                0.0090946163001294680, -0.0068367523674594463, -0.0013145571502701280,
                                                -0.00077550681812703609, -0.73615185577991682], dtype=np.float64)

CUTOFF_LATITUDE_DUSKSIDE_LEFT_INTERCEPT = np.float64(66.049010100315940)

CUTOFF_LATITUDE_DUSKSIDE_RIGHT_COEFFS = np.array([-0.036245424584295020, 0.017071209759156603, -0.85924889193050080,
                                                 0.016462770845528787, -0.053079177781167323, -0.0039780175521702786,
                                                 -0.0016200943643020558, -0.89500123807686616], dtype=np.float64)

CUTOFF_LATITUDE_DUSKSIDE_RIGHT_INTERCEPT = np.float64(70.673034286676597)

CUTOFF_LATITUDE_DUSKSIDE_BREAKPOINT = np.float64(-0.64624644044364477)

######################################################################

CUTOFF_LATITUDE_NIGHTSIDE_LEFT_COEFFS = np.array([-0.017353049858744189, 0.0064537525807491607, -0.71096166282243090,
                                                 0.0062946125421443387, -0.021552508655117458, -0.0015433462003364598,
                                                 -0.00074391041629426596, -0.52745013100999927], dtype=np.float64)

CUTOFF_LATITUDE_NIGHTSIDE_LEFT_INTERCEPT = np.float64(65.730947619280826)

CUTOFF_LATITUDE_NIGHTSIDE_RIGHT_COEFFS = np.array([-0.036156513215120617, 0.012544658847798351, -1.4175298304827511,
                                                  0.012199929144102633, -0.067383817784419098, -0.0042931701062772498,
                                                  -0.0015632046449290387, -0.61377937122507475], dtype=np.float64)

CUTOFF_LATITUDE_NIGHTSIDE_RIGHT_INTERCEPT = np.float64(70.130439080931808)

CUTOFF_LATITUDE_NIGHTSIDE_BREAKPOINT = np.float64(-0.28667911139357116)

#######################################################################

CUTOFF_LATITUDE_DAWNSIDE_LEFT_COEFFS = np.array([0.0026774439261010759, 0.0065105271857474421, 0.13523658681008649,
                                                0.0054509109717034406, 0.018763750923900906, -0.00011970967812424032,
                                                -0.00024716350157779401, -1.0651695579037004], dtype=np.float64)

CUTOFF_LATITUDE_DAWNSIDE_LEFT_INTERCEPT = np.float64(63.760884387214382)

CUTOFF_LATITUDE_DAWNSIDE_RIGHT_COEFFS = np.array([-0.012828416379424702, 0.014896359349719862, -0.36148878240305260,
                                                 0.013084585991682104, -0.019699253488552601, -0.0026607746283272143,
                                                 -0.0011258574203175681, -1.4269318707683063], dtype=np.float64)

CUTOFF_LATITUDE_DAWNSIDE_RIGHT_INTERCEPT = np.float64(71.376809209899861)

CUTOFF_LATITUDE_DAWNSIDE_BREAKPOINT = np.float64(-3.4002935935762295)

########################################################################

def _calculate_pls1_dayside(Kp, Dst, Pdyn, SYM_H, Bt, Vt, AE, proton_energy):
    return (PLS1_DAYSIDE_COEFFS[0] * Kp + PLS1_DAYSIDE_COEFFS[1] * Dst + PLS1_DAYSIDE_COEFFS[2] * Pdyn ** (1/3)
           + PLS1_DAYSIDE_COEFFS[3] * SYM_H + PLS1_DAYSIDE_COEFFS[4] * Bt + PLS1_DAYSIDE_COEFFS[5] * Vt
           + PLS1_DAYSIDE_COEFFS[6] * AE +  PLS1_DAYSIDE_COEFFS[7] * np.log(proton_energy) + PLS1_DAYSIDE_INTERCEPT)

def _calculate_pls1_duskside(Kp, Dst, Pdyn, SYM_H, Bt, Vt, AE, proton_energy):
    return (PLS1_DUSKSIDE_COEFFS[0] * Kp + PLS1_DUSKSIDE_COEFFS[1] * Dst + PLS1_DUSKSIDE_COEFFS[2] * Pdyn ** (1/3)
           + PLS1_DUSKSIDE_COEFFS[3] * SYM_H + PLS1_DUSKSIDE_COEFFS[4] * Bt + PLS1_DUSKSIDE_COEFFS[5] * Vt
           + PLS1_DUSKSIDE_COEFFS[6] * AE +  PLS1_DUSKSIDE_COEFFS[7] * np.log(proton_energy) + PLS1_DUSKSIDE_INTERCEPT)

def _calculate_pls1_nightside(Kp, Dst, Pdyn, SYM_H, Bt, Vt, AE, proton_energy):
    return (PLS1_NIGHTSIDE_COEFFS[0] * Kp + PLS1_NIGHTSIDE_COEFFS[1] * Dst + PLS1_NIGHTSIDE_COEFFS[2] * Pdyn ** (1/3)
           + PLS1_NIGHTSIDE_COEFFS[3] * SYM_H + PLS1_NIGHTSIDE_COEFFS[4] * Bt + PLS1_NIGHTSIDE_COEFFS[5] * Vt
           + PLS1_NIGHTSIDE_COEFFS[6] * AE +  PLS1_NIGHTSIDE_COEFFS[7] * np.log(proton_energy) + PLS1_NIGHTSIDE_INTERCEPT)

def _calculate_pls1_dawnside(Kp, Dst, Pdyn, SYM_H, Bt, Vt, AE, proton_energy):
    return (PLS1_DAWNSIDE_COEFFS[0] * Kp + PLS1_DAWNSIDE_COEFFS[1] * Dst + PLS1_DAWNSIDE_COEFFS[2] * Pdyn ** (1/3)
           + PLS1_DAWNSIDE_COEFFS[3] * SYM_H + PLS1_DAWNSIDE_COEFFS[4] * Bt + PLS1_DAWNSIDE_COEFFS[5] * Vt
           + PLS1_DAWNSIDE_COEFFS[6] * AE +  PLS1_DAWNSIDE_COEFFS[7] * np.log(proton_energy) + PLS1_DAWNSIDE_INTERCEPT)

##########################################################################

def calculate_cutoff_latitude_dayside(Kp, Dst, Pdyn, SYM_H, Bt, Vt, AE, proton_energy):
    pls1_dayside = _calculate_pls1_dayside(Kp, Dst, Pdyn, SYM_H, Bt, Vt, AE, proton_energy)
    
    cutoff_latitude_dayside = np.full_like(pls1_dayside, np.nan, dtype=np.float64)

    valid_mask = ~np.isnan(pls1_dayside)
    left_mask = valid_mask & (pls1_dayside < CUTOFF_LATITUDE_DAYSIDE_BREAKPOINT)
    right_mask = valid_mask & (pls1_dayside >= CUTOFF_LATITUDE_DAYSIDE_BREAKPOINT)
    
    cutoff_latitude_dayside[left_mask] = (
        CUTOFF_LATITUDE_DAYSIDE_LEFT_COEFFS[0] * Kp[left_mask] + CUTOFF_LATITUDE_DAYSIDE_LEFT_COEFFS[1] * Dst[left_mask]
        + CUTOFF_LATITUDE_DAYSIDE_LEFT_COEFFS[2] * (Pdyn ** (1/3))[left_mask] + CUTOFF_LATITUDE_DAYSIDE_LEFT_COEFFS[3] * SYM_H[left_mask]
        + CUTOFF_LATITUDE_DAYSIDE_LEFT_COEFFS[4] * Bt[left_mask] + CUTOFF_LATITUDE_DAYSIDE_LEFT_COEFFS[5] * Vt[left_mask]
        + CUTOFF_LATITUDE_DAYSIDE_LEFT_COEFFS[6] * AE[left_mask] + CUTOFF_LATITUDE_DAYSIDE_LEFT_COEFFS[7] * (np.log(proton_energy))[left_mask]
        + CUTOFF_LATITUDE_DAYSIDE_LEFT_INTERCEPT
    )
    
    cutoff_latitude_dayside[right_mask] = (
        CUTOFF_LATITUDE_DAYSIDE_RIGHT_COEFFS[0] * Kp[right_mask] + CUTOFF_LATITUDE_DAYSIDE_RIGHT_COEFFS[1] * Dst[right_mask]
        + CUTOFF_LATITUDE_DAYSIDE_RIGHT_COEFFS[2] * (Pdyn ** (1/3))[right_mask] + CUTOFF_LATITUDE_DAYSIDE_RIGHT_COEFFS[3] * SYM_H[right_mask]
        + CUTOFF_LATITUDE_DAYSIDE_RIGHT_COEFFS[4] * Bt[right_mask] + CUTOFF_LATITUDE_DAYSIDE_RIGHT_COEFFS[5] * Vt[right_mask]
        + CUTOFF_LATITUDE_DAYSIDE_RIGHT_COEFFS[6] * AE[right_mask] + CUTOFF_LATITUDE_DAYSIDE_RIGHT_COEFFS[7] * (np.log(proton_energy))[right_mask]
        + CUTOFF_LATITUDE_DAYSIDE_RIGHT_INTERCEPT
    )
    
    return cutoff_latitude_dayside

#############################################################################

def calculate_cutoff_latitude_duskside(Kp, Dst, Pdyn, SYM_H, Bt, Vt, AE, proton_energy):
    pls1_duskside = _calculate_pls1_duskside(Kp, Dst, Pdyn, SYM_H, Bt, Vt, AE, proton_energy)
    
    cutoff_latitude_duskside = np.full_like(pls1_duskside, np.nan, dtype=np.float64)

    valid_mask = ~np.isnan(pls1_duskside)
    left_mask = valid_mask & (pls1_duskside < CUTOFF_LATITUDE_DUSKSIDE_BREAKPOINT)
    right_mask = valid_mask & (pls1_duskside >= CUTOFF_LATITUDE_DUSKSIDE_BREAKPOINT)
    
    cutoff_latitude_duskside[left_mask] = (
        CUTOFF_LATITUDE_DUSKSIDE_LEFT_COEFFS[0] * Kp[left_mask] + CUTOFF_LATITUDE_DUSKSIDE_LEFT_COEFFS[1] * Dst[left_mask]
        + CUTOFF_LATITUDE_DUSKSIDE_LEFT_COEFFS[2] * (Pdyn ** (1/3))[left_mask] + CUTOFF_LATITUDE_DUSKSIDE_LEFT_COEFFS[3] * SYM_H[left_mask]
        + CUTOFF_LATITUDE_DUSKSIDE_LEFT_COEFFS[4] * Bt[left_mask] + CUTOFF_LATITUDE_DUSKSIDE_LEFT_COEFFS[5] * Vt[left_mask]
        + CUTOFF_LATITUDE_DUSKSIDE_LEFT_COEFFS[6] * AE[left_mask] + CUTOFF_LATITUDE_DUSKSIDE_LEFT_COEFFS[7] * (np.log(proton_energy))[left_mask]
        + CUTOFF_LATITUDE_DUSKSIDE_LEFT_INTERCEPT
    )
    
    cutoff_latitude_duskside[right_mask] = (
        CUTOFF_LATITUDE_DUSKSIDE_RIGHT_COEFFS[0] * Kp[right_mask] + CUTOFF_LATITUDE_DUSKSIDE_RIGHT_COEFFS[1] * Dst[right_mask]
        + CUTOFF_LATITUDE_DUSKSIDE_RIGHT_COEFFS[2] * (Pdyn ** (1/3))[right_mask] + CUTOFF_LATITUDE_DUSKSIDE_RIGHT_COEFFS[3] * SYM_H[right_mask]
        + CUTOFF_LATITUDE_DUSKSIDE_RIGHT_COEFFS[4] * Bt[right_mask] + CUTOFF_LATITUDE_DUSKSIDE_RIGHT_COEFFS[5] * Vt[right_mask]
        + CUTOFF_LATITUDE_DUSKSIDE_RIGHT_COEFFS[6] * AE[right_mask] + CUTOFF_LATITUDE_DUSKSIDE_RIGHT_COEFFS[7] * (np.log(proton_energy))[right_mask]
        + CUTOFF_LATITUDE_DUSKSIDE_RIGHT_INTERCEPT
    )
    
    return cutoff_latitude_duskside

#################################################################################

def calculate_cutoff_latitude_nightside(Kp, Dst, Pdyn, SYM_H, Bt, Vt, AE, proton_energy):
    pls1_nightside = _calculate_pls1_nightside(Kp, Dst, Pdyn, SYM_H, Bt, Vt, AE, proton_energy)
    
    cutoff_latitude_nightside = np.full_like(pls1_nightside, np.nan, dtype=np.float64)

    valid_mask = ~np.isnan(pls1_nightside)
    left_mask = valid_mask & (pls1_nightside < CUTOFF_LATITUDE_NIGHTSIDE_BREAKPOINT)
    right_mask = valid_mask & (pls1_nightside >= CUTOFF_LATITUDE_NIGHTSIDE_BREAKPOINT)
    
    cutoff_latitude_nightside[left_mask] = (
        CUTOFF_LATITUDE_NIGHTSIDE_LEFT_COEFFS[0] * Kp[left_mask] + CUTOFF_LATITUDE_NIGHTSIDE_LEFT_COEFFS[1] * Dst[left_mask]
        + CUTOFF_LATITUDE_NIGHTSIDE_LEFT_COEFFS[2] * (Pdyn ** (1/3))[left_mask] + CUTOFF_LATITUDE_NIGHTSIDE_LEFT_COEFFS[3] * SYM_H[left_mask]
        + CUTOFF_LATITUDE_NIGHTSIDE_LEFT_COEFFS[4] * Bt[left_mask] + CUTOFF_LATITUDE_NIGHTSIDE_LEFT_COEFFS[5] * Vt[left_mask]
        + CUTOFF_LATITUDE_NIGHTSIDE_LEFT_COEFFS[6] * AE[left_mask] + CUTOFF_LATITUDE_NIGHTSIDE_LEFT_COEFFS[7] * (np.log(proton_energy))[left_mask]
        + CUTOFF_LATITUDE_NIGHTSIDE_LEFT_INTERCEPT
    )
    
    cutoff_latitude_nightside[right_mask] = (
        CUTOFF_LATITUDE_NIGHTSIDE_RIGHT_COEFFS[0] * Kp[right_mask] + CUTOFF_LATITUDE_NIGHTSIDE_RIGHT_COEFFS[1] * Dst[right_mask]
        + CUTOFF_LATITUDE_NIGHTSIDE_RIGHT_COEFFS[2] * (Pdyn ** (1/3))[right_mask] + CUTOFF_LATITUDE_NIGHTSIDE_RIGHT_COEFFS[3] * SYM_H[right_mask]
        + CUTOFF_LATITUDE_NIGHTSIDE_RIGHT_COEFFS[4] * Bt[right_mask] + CUTOFF_LATITUDE_NIGHTSIDE_RIGHT_COEFFS[5] * Vt[right_mask]
        + CUTOFF_LATITUDE_NIGHTSIDE_RIGHT_COEFFS[6] * AE[right_mask] + CUTOFF_LATITUDE_NIGHTSIDE_RIGHT_COEFFS[7] * (np.log(proton_energy))[right_mask]
        + CUTOFF_LATITUDE_NIGHTSIDE_RIGHT_INTERCEPT
    )
    
    return cutoff_latitude_nightside

##################################################################################

def calculate_cutoff_latitude_dawnside(Kp, Dst, Pdyn, SYM_H, Bt, Vt, AE, proton_energy):
    pls1_dawnside = _calculate_pls1_dawnside(Kp, Dst, Pdyn, SYM_H, Bt, Vt, AE, proton_energy)
    
    cutoff_latitude_dawnside = np.full_like(pls1_dawnside, np.nan, dtype=np.float64)

    valid_mask = ~np.isnan(pls1_dawnside)
    left_mask = valid_mask & (pls1_dawnside < CUTOFF_LATITUDE_DAWNSIDE_BREAKPOINT)
    right_mask = valid_mask & (pls1_dawnside >= CUTOFF_LATITUDE_DAWNSIDE_BREAKPOINT)
    
    cutoff_latitude_dawnside[left_mask] = (
        CUTOFF_LATITUDE_DAWNSIDE_LEFT_COEFFS[0] * Kp[left_mask] + CUTOFF_LATITUDE_DAWNSIDE_LEFT_COEFFS[1] * Dst[left_mask]
        + CUTOFF_LATITUDE_DAWNSIDE_LEFT_COEFFS[2] * (Pdyn ** (1/3))[left_mask] + CUTOFF_LATITUDE_DAWNSIDE_LEFT_COEFFS[3] * SYM_H[left_mask]
        + CUTOFF_LATITUDE_DAWNSIDE_LEFT_COEFFS[4] * Bt[left_mask] + CUTOFF_LATITUDE_DAWNSIDE_LEFT_COEFFS[5] * Vt[left_mask]
        + CUTOFF_LATITUDE_DAWNSIDE_LEFT_COEFFS[6] * AE[left_mask] + CUTOFF_LATITUDE_DAWNSIDE_LEFT_COEFFS[7] * (np.log(proton_energy))[left_mask]
        + CUTOFF_LATITUDE_DAWNSIDE_LEFT_INTERCEPT
    )
    
    cutoff_latitude_dawnside[right_mask] = (
        CUTOFF_LATITUDE_DAWNSIDE_RIGHT_COEFFS[0] * Kp[right_mask] + CUTOFF_LATITUDE_DAWNSIDE_RIGHT_COEFFS[1] * Dst[right_mask]
        + CUTOFF_LATITUDE_DAWNSIDE_RIGHT_COEFFS[2] * (Pdyn ** (1/3))[right_mask] + CUTOFF_LATITUDE_DAWNSIDE_RIGHT_COEFFS[3] * SYM_H[right_mask]
        + CUTOFF_LATITUDE_DAWNSIDE_RIGHT_COEFFS[4] * Bt[right_mask] + CUTOFF_LATITUDE_DAWNSIDE_RIGHT_COEFFS[5] * Vt[right_mask]
        + CUTOFF_LATITUDE_DAWNSIDE_RIGHT_COEFFS[6] * AE[right_mask] + CUTOFF_LATITUDE_DAWNSIDE_RIGHT_COEFFS[7] * (np.log(proton_energy))[right_mask]
        + CUTOFF_LATITUDE_DAWNSIDE_RIGHT_INTERCEPT
    )
    
    return cutoff_latitude_dawnside

####################################################################################

def calculate_cutoff_latitude_divisional_mlt(Kp, Dst, Pdyn, SYM_H, Bt, Vt, AE, proton_energy, MLT):
    cutoff_latitude_divisional_mlt = np.full_like(Kp, np.nan, dtype=np.float64)

    dayside_mask = (MLT >= 9.0) & (MLT < 15.0)

    duskside_mask = (MLT >= 15.0) & (MLT < 21.0)

    nightside_mask = ((MLT >= 21.0) & (MLT <= 24.0)) | ((MLT >= 0.0) & (MLT < 3.0))

    dawnside_mask = (MLT >= 3.0) & (MLT < 9.0)

    if np.any(dayside_mask):
        cutoff_latitude_divisional_mlt[dayside_mask] = calculate_cutoff_latitude_dayside(
            Kp[dayside_mask],
            Dst[dayside_mask],
            Pdyn[dayside_mask],
            SYM_H[dayside_mask],
            Bt[dayside_mask],
            Vt[dayside_mask],
            AE[dayside_mask],
            proton_energy[dayside_mask],
        )

    if np.any(duskside_mask):
        cutoff_latitude_divisional_mlt[duskside_mask] = calculate_cutoff_latitude_duskside(
            Kp[duskside_mask],
            Dst[duskside_mask],
            Pdyn[duskside_mask],
            SYM_H[duskside_mask],
            Bt[duskside_mask],
            Vt[duskside_mask],
            AE[duskside_mask],
            proton_energy[duskside_mask],
        )

    if np.any(nightside_mask):
        cutoff_latitude_divisional_mlt[nightside_mask] = calculate_cutoff_latitude_nightside(
            Kp[nightside_mask],
            Dst[nightside_mask],
            Pdyn[nightside_mask],
            SYM_H[nightside_mask],
            Bt[nightside_mask],
            Vt[nightside_mask],
            AE[nightside_mask],
            proton_energy[nightside_mask],
        )

    if np.any(dawnside_mask):
        cutoff_latitude_divisional_mlt[dawnside_mask] = calculate_cutoff_latitude_dawnside(
            Kp[dawnside_mask],
            Dst[dawnside_mask],
            Pdyn[dawnside_mask],
            SYM_H[dawnside_mask],
            Bt[dawnside_mask],
            Vt[dawnside_mask],
            AE[dawnside_mask],
            proton_energy[dawnside_mask],
        )

    return cutoff_latitude_divisional_mlt

##############################################################################

def calculate_cutoff_latitude_continuous_mlt(Kp, Dst, Pdyn, SYM_H, Bt, Vt, AE, proton_energy, MLT):
    cutoff_latitude_12mlt = calculate_cutoff_latitude_dayside(
        Kp, Dst, Pdyn, SYM_H, Bt, Vt, AE, proton_energy
    )

    cutoff_latitude_18mlt = calculate_cutoff_latitude_duskside(
        Kp, Dst, Pdyn, SYM_H, Bt, Vt, AE, proton_energy
    )

    cutoff_latitude_00mlt = calculate_cutoff_latitude_nightside(
        Kp, Dst, Pdyn, SYM_H, Bt, Vt, AE, proton_energy
    )

    cutoff_latitude_06mlt = calculate_cutoff_latitude_dawnside(
        Kp, Dst, Pdyn, SYM_H, Bt, Vt, AE, proton_energy
    )

    cutoff_latitude_continuous_mlt = np.full_like(Kp, np.nan, dtype=np.float64)

    mask_00_to_06 = (MLT >= 0.0) & (MLT < 6.0)
    mask_06_to_12 = (MLT >= 6.0) & (MLT < 12.0)
    mask_12_to_18 = (MLT >= 12.0) & (MLT < 18.0)
    mask_18_to_24 = (MLT >= 18.0) & (MLT <= 24.0)

    cutoff_latitude_continuous_mlt[mask_00_to_06] = (
        cutoff_latitude_00mlt[mask_00_to_06]
        + (
            cutoff_latitude_06mlt[mask_00_to_06]
            - cutoff_latitude_00mlt[mask_00_to_06]
        )
        * MLT[mask_00_to_06] / 6.0
    )

    cutoff_latitude_continuous_mlt[mask_06_to_12] = (
        cutoff_latitude_06mlt[mask_06_to_12]
        + (
            cutoff_latitude_12mlt[mask_06_to_12]
            - cutoff_latitude_06mlt[mask_06_to_12]
        )
        * (MLT[mask_06_to_12] - 6.0) / 6.0
    )

    cutoff_latitude_continuous_mlt[mask_12_to_18] = (
        cutoff_latitude_12mlt[mask_12_to_18]
        + (
            cutoff_latitude_18mlt[mask_12_to_18]
            - cutoff_latitude_12mlt[mask_12_to_18]
        )
        * (MLT[mask_12_to_18] - 12.0) / 6.0
    )

    cutoff_latitude_continuous_mlt[mask_18_to_24] = (
        cutoff_latitude_18mlt[mask_18_to_24]
        + (
            cutoff_latitude_00mlt[mask_18_to_24]
            - cutoff_latitude_18mlt[mask_18_to_24]
        )
        * (MLT[mask_18_to_24] - 18.0) / 6.0
    )

    return cutoff_latitude_continuous_mlt

#######################################################################

def _calculate_pls1_non_energy_term(
    Kp, Dst, Pdyn, SYM_H, Bt, Vt, AE, pls1_coeffs, pls1_intercept
):
    return (
        pls1_coeffs[0] * Kp
        + pls1_coeffs[1] * Dst
        + pls1_coeffs[2] * Pdyn ** (1.0 / 3.0)
        + pls1_coeffs[3] * SYM_H
        + pls1_coeffs[4] * Bt
        + pls1_coeffs[5] * Vt
        + pls1_coeffs[6] * AE
        + pls1_intercept
    )

#######################################################################

def _calculate_cutoff_latitude_non_energy_term(
    Kp, Dst, Pdyn, SYM_H, Bt, Vt, AE,
    cutoff_latitude_coeffs, cutoff_latitude_intercept
):
    return (
        cutoff_latitude_coeffs[0] * Kp
        + cutoff_latitude_coeffs[1] * Dst
        + cutoff_latitude_coeffs[2] * Pdyn ** (1.0 / 3.0)
        + cutoff_latitude_coeffs[3] * SYM_H
        + cutoff_latitude_coeffs[4] * Bt
        + cutoff_latitude_coeffs[5] * Vt
        + cutoff_latitude_coeffs[6] * AE
        + cutoff_latitude_intercept
    )

#######################################################################

def _calculate_proton_cutoff_energy_single_sector(
    Kp, Dst, Pdyn, SYM_H, Bt, Vt, AE, Latitude,
    pls1_coeffs, pls1_intercept,
    left_coeffs, left_intercept,
    right_coeffs, right_intercept,
    breakpoint,
):
    log_proton_cutoff_energy = np.full_like(Kp, np.nan, dtype=np.float64)

    pls1_non_energy = _calculate_pls1_non_energy_term(
        Kp, Dst, Pdyn, SYM_H, Bt, Vt, AE,
        pls1_coeffs, pls1_intercept
    )

    pls1_log_energy_coeff = pls1_coeffs[7]

    left_non_energy = _calculate_cutoff_latitude_non_energy_term(
        Kp, Dst, Pdyn, SYM_H, Bt, Vt, AE,
        left_coeffs, left_intercept
    )

    right_non_energy = _calculate_cutoff_latitude_non_energy_term(
        Kp, Dst, Pdyn, SYM_H, Bt, Vt, AE,
        right_coeffs, right_intercept
    )

    left_log_energy_coeff = left_coeffs[7]
    right_log_energy_coeff = right_coeffs[7]

    log_energy_left = (Latitude - left_non_energy) / left_log_energy_coeff

    log_energy_right = (Latitude - right_non_energy) / right_log_energy_coeff

    pls1_left = pls1_non_energy + pls1_log_energy_coeff * log_energy_left
    pls1_right = pls1_non_energy + pls1_log_energy_coeff * log_energy_right

    valid_left = np.isfinite(log_energy_left) & (pls1_left < breakpoint)
    valid_right = np.isfinite(log_energy_right) & (pls1_right >= breakpoint)

    log_proton_cutoff_energy[valid_left] = log_energy_left[valid_left]
    log_proton_cutoff_energy[valid_right] = log_energy_right[valid_right]
    
    proton_cutoff_energy = np.exp(log_proton_cutoff_energy)

    return proton_cutoff_energy

######################################################################

def calculate_proton_cutoff_energy_dayside(Kp, Dst, Pdyn, SYM_H, Bt, Vt, AE, Latitude):
    return _calculate_proton_cutoff_energy_single_sector(
        Kp, Dst, Pdyn, SYM_H, Bt, Vt, AE, Latitude,
        PLS1_DAYSIDE_COEFFS,
        PLS1_DAYSIDE_INTERCEPT,
        CUTOFF_LATITUDE_DAYSIDE_LEFT_COEFFS,
        CUTOFF_LATITUDE_DAYSIDE_LEFT_INTERCEPT,
        CUTOFF_LATITUDE_DAYSIDE_RIGHT_COEFFS,
        CUTOFF_LATITUDE_DAYSIDE_RIGHT_INTERCEPT,
        CUTOFF_LATITUDE_DAYSIDE_BREAKPOINT
    )

#######################################################################

def calculate_proton_cutoff_energy_duskside(Kp, Dst, Pdyn, SYM_H, Bt, Vt, AE, Latitude):
    return _calculate_proton_cutoff_energy_single_sector(
        Kp, Dst, Pdyn, SYM_H, Bt, Vt, AE, Latitude,
        PLS1_DUSKSIDE_COEFFS,
        PLS1_DUSKSIDE_INTERCEPT,
        CUTOFF_LATITUDE_DUSKSIDE_LEFT_COEFFS,
        CUTOFF_LATITUDE_DUSKSIDE_LEFT_INTERCEPT,
        CUTOFF_LATITUDE_DUSKSIDE_RIGHT_COEFFS,
        CUTOFF_LATITUDE_DUSKSIDE_RIGHT_INTERCEPT,
        CUTOFF_LATITUDE_DUSKSIDE_BREAKPOINT
    )

#######################################################################

def calculate_proton_cutoff_energy_nightside(Kp, Dst, Pdyn, SYM_H, Bt, Vt, AE, Latitude):
    return _calculate_proton_cutoff_energy_single_sector(
        Kp, Dst, Pdyn, SYM_H, Bt, Vt, AE, Latitude,
        PLS1_NIGHTSIDE_COEFFS,
        PLS1_NIGHTSIDE_INTERCEPT,
        CUTOFF_LATITUDE_NIGHTSIDE_LEFT_COEFFS,
        CUTOFF_LATITUDE_NIGHTSIDE_LEFT_INTERCEPT,
        CUTOFF_LATITUDE_NIGHTSIDE_RIGHT_COEFFS,
        CUTOFF_LATITUDE_NIGHTSIDE_RIGHT_INTERCEPT,
        CUTOFF_LATITUDE_NIGHTSIDE_BREAKPOINT
    )

########################################################################

def calculate_proton_cutoff_energy_dawnside(Kp, Dst, Pdyn, SYM_H, Bt, Vt, AE, Latitude):
    return _calculate_proton_cutoff_energy_single_sector(
        Kp, Dst, Pdyn, SYM_H, Bt, Vt, AE, Latitude,
        PLS1_DAWNSIDE_COEFFS,
        PLS1_DAWNSIDE_INTERCEPT,
        CUTOFF_LATITUDE_DAWNSIDE_LEFT_COEFFS,
        CUTOFF_LATITUDE_DAWNSIDE_LEFT_INTERCEPT,
        CUTOFF_LATITUDE_DAWNSIDE_RIGHT_COEFFS,
        CUTOFF_LATITUDE_DAWNSIDE_RIGHT_INTERCEPT,
        CUTOFF_LATITUDE_DAWNSIDE_BREAKPOINT
    )

#########################################################################

def calculate_proton_cutoff_energy_divisional_mlt(
    Kp, Dst, Pdyn, SYM_H, Bt, Vt, AE, Latitude, MLT
):
    proton_cutoff_energy_divisional_mlt = np.full_like(Kp, np.nan, dtype=np.float64)

    dayside_mask = (MLT >= 9.0) & (MLT < 15.0)
    duskside_mask = (MLT >= 15.0) & (MLT < 21.0)
    nightside_mask = ((MLT >= 21.0) & (MLT <= 24.0)) | ((MLT >= 0.0) & (MLT < 3.0))
    dawnside_mask = (MLT >= 3.0) & (MLT < 9.0)

    if np.any(dayside_mask):
        proton_cutoff_energy_divisional_mlt[dayside_mask] = calculate_proton_cutoff_energy_dayside(
            Kp[dayside_mask],
            Dst[dayside_mask],
            Pdyn[dayside_mask],
            SYM_H[dayside_mask],
            Bt[dayside_mask],
            Vt[dayside_mask],
            AE[dayside_mask],
            Latitude[dayside_mask],
        )

    if np.any(duskside_mask):
        proton_cutoff_energy_divisional_mlt[duskside_mask] = calculate_proton_cutoff_energy_duskside(
            Kp[duskside_mask],
            Dst[duskside_mask],
            Pdyn[duskside_mask],
            SYM_H[duskside_mask],
            Bt[duskside_mask],
            Vt[duskside_mask],
            AE[duskside_mask],
            Latitude[duskside_mask],
        )

    if np.any(nightside_mask):
        proton_cutoff_energy_divisional_mlt[nightside_mask] = calculate_proton_cutoff_energy_nightside(
            Kp[nightside_mask],
            Dst[nightside_mask],
            Pdyn[nightside_mask],
            SYM_H[nightside_mask],
            Bt[nightside_mask],
            Vt[nightside_mask],
            AE[nightside_mask],
            Latitude[nightside_mask],
        )

    if np.any(dawnside_mask):
        proton_cutoff_energy_divisional_mlt[dawnside_mask] = calculate_proton_cutoff_energy_dawnside(
            Kp[dawnside_mask],
            Dst[dawnside_mask],
            Pdyn[dawnside_mask],
            SYM_H[dawnside_mask],
            Bt[dawnside_mask],
            Vt[dawnside_mask],
            AE[dawnside_mask],
            Latitude[dawnside_mask],
        )

    return proton_cutoff_energy_divisional_mlt

######################################################################

def _get_sector_parameters(sector_name):
    if sector_name == "dayside":
        return (
            PLS1_DAYSIDE_COEFFS,
            PLS1_DAYSIDE_INTERCEPT,
            CUTOFF_LATITUDE_DAYSIDE_LEFT_COEFFS,
            CUTOFF_LATITUDE_DAYSIDE_LEFT_INTERCEPT,
            CUTOFF_LATITUDE_DAYSIDE_RIGHT_COEFFS,
            CUTOFF_LATITUDE_DAYSIDE_RIGHT_INTERCEPT,
            CUTOFF_LATITUDE_DAYSIDE_BREAKPOINT,
        )

    if sector_name == "duskside":
        return (
            PLS1_DUSKSIDE_COEFFS,
            PLS1_DUSKSIDE_INTERCEPT,
            CUTOFF_LATITUDE_DUSKSIDE_LEFT_COEFFS,
            CUTOFF_LATITUDE_DUSKSIDE_LEFT_INTERCEPT,
            CUTOFF_LATITUDE_DUSKSIDE_RIGHT_COEFFS,
            CUTOFF_LATITUDE_DUSKSIDE_RIGHT_INTERCEPT,
            CUTOFF_LATITUDE_DUSKSIDE_BREAKPOINT,
        )

    if sector_name == "nightside":
        return (
            PLS1_NIGHTSIDE_COEFFS,
            PLS1_NIGHTSIDE_INTERCEPT,
            CUTOFF_LATITUDE_NIGHTSIDE_LEFT_COEFFS,
            CUTOFF_LATITUDE_NIGHTSIDE_LEFT_INTERCEPT,
            CUTOFF_LATITUDE_NIGHTSIDE_RIGHT_COEFFS,
            CUTOFF_LATITUDE_NIGHTSIDE_RIGHT_INTERCEPT,
            CUTOFF_LATITUDE_NIGHTSIDE_BREAKPOINT,
        )

    if sector_name == "dawnside":
        return (
            PLS1_DAWNSIDE_COEFFS,
            PLS1_DAWNSIDE_INTERCEPT,
            CUTOFF_LATITUDE_DAWNSIDE_LEFT_COEFFS,
            CUTOFF_LATITUDE_DAWNSIDE_LEFT_INTERCEPT,
            CUTOFF_LATITUDE_DAWNSIDE_RIGHT_COEFFS,
            CUTOFF_LATITUDE_DAWNSIDE_RIGHT_INTERCEPT,
            CUTOFF_LATITUDE_DAWNSIDE_BREAKPOINT,
        )

    raise ValueError("Unknown sector_name.")

####################################################################

def _calculate_sector_parts(Kp, Dst, Pdyn, SYM_H, Bt, Vt, AE, sector_name):
    (
        pls1_coeffs,
        pls1_intercept,
        left_coeffs,
        left_intercept,
        right_coeffs,
        right_intercept,
        breakpoint,
    ) = _get_sector_parameters(sector_name)

    pls1_non_energy = _calculate_pls1_non_energy_term(
        Kp, Dst, Pdyn, SYM_H, Bt, Vt, AE,
        pls1_coeffs, pls1_intercept
    )

    left_non_energy = _calculate_cutoff_latitude_non_energy_term(
        Kp, Dst, Pdyn, SYM_H, Bt, Vt, AE,
        left_coeffs, left_intercept
    )

    right_non_energy = _calculate_cutoff_latitude_non_energy_term(
        Kp, Dst, Pdyn, SYM_H, Bt, Vt, AE,
        right_coeffs, right_intercept
    )

    return {
        "pls1_non_energy": pls1_non_energy,
        "pls1_log_energy_coeff": pls1_coeffs[7],
        "left_non_energy": left_non_energy,
        "left_log_energy_coeff": left_coeffs[7],
        "right_non_energy": right_non_energy,
        "right_log_energy_coeff": right_coeffs[7],
        "breakpoint": breakpoint,
    }

#####################################################################

def _solve_proton_cutoff_energy_continuous_between_two_sectors(
    Kp, Dst, Pdyn, SYM_H, Bt, Vt, AE, Latitude,
    weight_first_sector, weight_second_sector,
    first_sector_name, second_sector_name
):
    log_proton_cutoff_energy = np.full_like(Kp, np.nan, dtype=np.float64)

    first = _calculate_sector_parts(
        Kp, Dst, Pdyn, SYM_H, Bt, Vt, AE, first_sector_name
    )

    second = _calculate_sector_parts(
        Kp, Dst, Pdyn, SYM_H, Bt, Vt, AE, second_sector_name
    )

    sector_branches = ("left", "right")

    first_active = weight_first_sector > 0.0
    second_active = weight_second_sector > 0.0

    for first_branch in sector_branches:
        for second_branch in sector_branches:
            first_non_energy = first[f"{first_branch}_non_energy"]
            first_log_energy_coeff = first[f"{first_branch}_log_energy_coeff"]

            second_non_energy = second[f"{second_branch}_non_energy"]
            second_log_energy_coeff = second[f"{second_branch}_log_energy_coeff"]

            combined_non_energy = (
                weight_first_sector * first_non_energy
                + weight_second_sector * second_non_energy
            )

            combined_log_energy_coeff = (
                weight_first_sector * first_log_energy_coeff
                + weight_second_sector * second_log_energy_coeff
            )

            log_energy_candidate = (
                Latitude
                - combined_non_energy
            ) / combined_log_energy_coeff

            first_pls1_candidate = (
                first["pls1_non_energy"]
                + first["pls1_log_energy_coeff"] * log_energy_candidate
            )

            second_pls1_candidate = (
                second["pls1_non_energy"]
                + second["pls1_log_energy_coeff"] * log_energy_candidate
            )

            if first_branch == "left":
                first_branch_valid = first_pls1_candidate < first["breakpoint"]
            else:
                first_branch_valid = first_pls1_candidate >= first["breakpoint"]

            if second_branch == "left":
                second_branch_valid = second_pls1_candidate < second["breakpoint"]
            else:
                second_branch_valid = second_pls1_candidate >= second["breakpoint"]

            # If one sector has exactly zero weight, its branch condition should not
            # reject the solution because that sector does not contribute to the
            # interpolated cutoff latitude.
            first_branch_valid = (~first_active) | first_branch_valid
            second_branch_valid = (~second_active) | second_branch_valid

            candidate_valid = (
                np.isfinite(log_energy_candidate)
                & first_branch_valid
                & second_branch_valid
            )

            fill_mask = candidate_valid & np.isnan(log_proton_cutoff_energy)

            log_proton_cutoff_energy[fill_mask] = log_energy_candidate[fill_mask]
            
            proton_cutoff_energy = np.exp(log_proton_cutoff_energy)

    return proton_cutoff_energy

#####################################################################

def calculate_proton_cutoff_energy_continuous_mlt(
    Kp, Dst, Pdyn, SYM_H, Bt, Vt, AE, Latitude, MLT
):
    proton_cutoff_energy_continuous_mlt = np.full_like(Kp, np.nan, dtype=np.float64)

    mask_00_to_06 = (MLT >= 0.0) & (MLT < 6.0)
    mask_06_to_12 = (MLT >= 6.0) & (MLT < 12.0)
    mask_12_to_18 = (MLT >= 12.0) & (MLT < 18.0)
    mask_18_to_24 = (MLT >= 18.0) & (MLT <= 24.0)

    if np.any(mask_00_to_06):
        weight_second_sector = MLT[mask_00_to_06] / 6.0
        weight_first_sector = 1.0 - weight_second_sector

        proton_energy = _solve_proton_cutoff_energy_continuous_between_two_sectors(
            Kp[mask_00_to_06],
            Dst[mask_00_to_06],
            Pdyn[mask_00_to_06],
            SYM_H[mask_00_to_06],
            Bt[mask_00_to_06],
            Vt[mask_00_to_06],
            AE[mask_00_to_06],
            Latitude[mask_00_to_06],
            weight_first_sector,
            weight_second_sector,
            "nightside",
            "dawnside",
        )

        proton_cutoff_energy_continuous_mlt[mask_00_to_06] = proton_energy

    if np.any(mask_06_to_12):
        weight_second_sector = (MLT[mask_06_to_12] - 6.0) / 6.0
        weight_first_sector = 1.0 - weight_second_sector

        proton_energy = _solve_proton_cutoff_energy_continuous_between_two_sectors(
            Kp[mask_06_to_12],
            Dst[mask_06_to_12],
            Pdyn[mask_06_to_12],
            SYM_H[mask_06_to_12],
            Bt[mask_06_to_12],
            Vt[mask_06_to_12],
            AE[mask_06_to_12],
            Latitude[mask_06_to_12],
            weight_first_sector,
            weight_second_sector,
            "dawnside",
            "dayside",
        )

        proton_cutoff_energy_continuous_mlt[mask_06_to_12] = proton_energy

    if np.any(mask_12_to_18):
        weight_second_sector = (MLT[mask_12_to_18] - 12.0) / 6.0
        weight_first_sector = 1.0 - weight_second_sector

        proton_energy = _solve_proton_cutoff_energy_continuous_between_two_sectors(
            Kp[mask_12_to_18],
            Dst[mask_12_to_18],
            Pdyn[mask_12_to_18],
            SYM_H[mask_12_to_18],
            Bt[mask_12_to_18],
            Vt[mask_12_to_18],
            AE[mask_12_to_18],
            Latitude[mask_12_to_18],
            weight_first_sector,
            weight_second_sector,
            "dayside",
            "duskside",
        )

        proton_cutoff_energy_continuous_mlt[mask_12_to_18] = proton_energy

    if np.any(mask_18_to_24):
        weight_second_sector = (MLT[mask_18_to_24] - 18.0) / 6.0
        weight_first_sector = 1.0 - weight_second_sector

        proton_energy = _solve_proton_cutoff_energy_continuous_between_two_sectors(
            Kp[mask_18_to_24],
            Dst[mask_18_to_24],
            Pdyn[mask_18_to_24],
            SYM_H[mask_18_to_24],
            Bt[mask_18_to_24],
            Vt[mask_18_to_24],
            AE[mask_18_to_24],
            Latitude[mask_18_to_24],
            weight_first_sector,
            weight_second_sector,
            "duskside",
            "nightside",
        )

        proton_cutoff_energy_continuous_mlt[mask_18_to_24] = proton_energy

    return proton_cutoff_energy_continuous_mlt
