# PLS-based Model of Cutoffs

This repository provides Python code for a PLS-based empirical model of cutoff latitude and proton cutoff energy for solar energetic protons.

The model is associated with the manuscript:

**A New Model of Cutoff Latitude for Solar Energetic Protons Based on the Partial Least Squares Method**

Authors:

Zhenghao She, Hui Zhu, Qian Song, Yao Chen, Yingying Zhao, Huicong Chen, and Zhijie Qin

## Description

This repository contains the Python implementation of a PLS-based model of cutoffs. The model can be used to calculate:

1. cutoff latitude from solar wind parameters, geomagnetic activity indices, proton energy, and magnetic local time;
2. proton cutoff energy by inverting the cutoff latitude formulation.

The main Python file is:

```text
PLS_based_model.py
```

The four main user-facing functions are:

```python
calculate_cutoff_latitude_divisional_mlt
calculate_cutoff_latitude_continuous_mlt
calculate_proton_cutoff_energy_divisional_mlt
calculate_proton_cutoff_energy_continuous_mlt
```

## Requirements

The code requires Python and NumPy.

```bash
pip install numpy
```

The only required Python package is:

```text
numpy
```

## Input variables for cutoff latitude calculation

The following two functions calculate cutoff latitude:

```python
calculate_cutoff_latitude_divisional_mlt(
    Kp, Dst, Pdyn, SYM_H, Bt, Vt, AE, proton_energy, MLT
)

calculate_cutoff_latitude_continuous_mlt(
    Kp, Dst, Pdyn, SYM_H, Bt, Vt, AE, proton_energy, MLT
)
```

All inputs must be one-dimensional NumPy arrays of floating-point values with the same length. Elements with the same index in all input arrays are assumed to correspond to the same time point.

| Variable        | Description                                        | Unit                 |
| --------------- | -------------------------------------------------- | -------------------- |
| `Kp`            | Kp index, in decimal representation, i.e., Kp × 10 | dimensionless        |
| `Dst`           | Dst index                                          | nT                   |
| `Pdyn`          | Solar wind dynamic pressure                        | nPa                  |
| `SYM_H`         | SYM-H index                                        | nT                   |
| `Bt`            | Total interplanetary magnetic field strength       | nT                   |
| `Vt`            | Solar wind speed                                   | km/s                 |
| `AE`            | AE index                                           | nT                   |
| `proton_energy` | Proton energy                                      | MeV                  |
| `MLT`           | Magnetic local time                                | hour, 0 <= MLT <= 24 |

The value of `proton_energy` must be positive because the model uses the natural logarithm of proton energy, i.e., `np.log(proton_energy)`.

The outputs are one-dimensional NumPy arrays of floating-point values:

```python
cutoff_latitude_divisional_mlt
cutoff_latitude_continuous_mlt
```

The unit of cutoff latitude is degree. The output arrays have the same length as the input arrays, and each output element corresponds to the same time point as the input elements with the same index.

## Input variables for proton cutoff energy calculation

The following two functions calculate proton cutoff energy:

```python
calculate_proton_cutoff_energy_divisional_mlt(
    Kp, Dst, Pdyn, SYM_H, Bt, Vt, AE, Latitude, MLT
)

calculate_proton_cutoff_energy_continuous_mlt(
    Kp, Dst, Pdyn, SYM_H, Bt, Vt, AE, Latitude, MLT
)
```

All inputs must be one-dimensional NumPy arrays of floating-point values with the same length. Elements with the same index in all input arrays are assumed to correspond to the same time point.

| Variable   | Description                                        | Unit                 |
| ---------- | -------------------------------------------------- | -------------------- |
| `Kp`       | Kp index, in decimal representation, i.e., Kp × 10 | dimensionless        |
| `Dst`      | Dst index                                          | nT                   |
| `Pdyn`     | Solar wind dynamic pressure                        | nPa                  |
| `SYM_H`    | SYM-H index                                        | nT                   |
| `Bt`       | Total interplanetary magnetic field strength       | nT                   |
| `Vt`       | Solar wind speed                                   | km/s                 |
| `AE`       | AE index                                           | nT                   |
| `Latitude` | Geomagnetic latitude                               | degree               |
| `MLT`      | Magnetic local time                                | hour, 0 <= MLT <= 24 |

The outputs are one-dimensional NumPy arrays of floating-point values:

```python
proton_cutoff_energy_divisional_mlt
proton_cutoff_energy_continuous_mlt
```

The unit of proton cutoff energy is MeV. The output arrays have the same length as the input arrays, and each output element corresponds to the same time point as the input elements with the same index.

## MLT treatment

The divisional MLT model uses four MLT sectors:

| Sector    | MLT range                        |
| --------- | -------------------------------- |
| dayside   | 9 <= MLT < 15                    |
| duskside  | 15 <= MLT < 21                   |
| nightside | 21 <= MLT <= 24 and 0 <= MLT < 3 |
| dawnside  | 3 <= MLT < 9                     |

The continuous MLT model treats the sectoral cutoff latitudes as values at four representative MLT locations:

| Representative MLT | Sectoral model |
| ------------------ | -------------- |
| 06 MLT             | dawnside       |
| 12 MLT             | dayside        |
| 18 MLT             | duskside       |
| 00/24 MLT          | nightside      |

The cutoff latitude is then obtained by linear interpolation in MLT between these representative locations. The continuous proton cutoff energy is obtained by inverting this continuous cutoff latitude formulation.

## Notes

All calculations are performed element by element. Therefore, all input arrays must have the same length and must be aligned in time.

Invalid or out-of-range values may produce `np.nan` in the output. In particular, MLT values outside the range `0 <= MLT <= 24` are not assigned to any MLT sector.

For cutoff latitude calculations, `proton_energy` must be positive because the model uses `np.log(proton_energy)`.

## Citation

If you use this model in research work intended for publication, please cite the manuscript associated with this model:

**A New Model of Cutoff Latitude for Solar Energetic Protons Based on the Partial Least Squares Method**

The final bibliographic information will be updated after the manuscript is published.

## Contact

For questions, bug reports, or problems encountered when using the model, please contact:

* Hui Zhu: [huizhu@email.sdu.edu.cn](mailto:huizhu@email.sdu.edu.cn)
* Zhenghao She: [shezhenghao@qq.com](mailto:shezhenghao@qq.com)

## Copyright

Copyright (c) 2026 Zhenghao She, Hui Zhu, Qian Song, Yao Chen, Yingying Zhao, Huicong Chen, and Zhijie Qin. All rights reserved.

## Repository link

The GitHub repository for this model is:

```text
https://github.com/zhenghaoshe/PLS-based-model-of-cutoffs
```

If the repository is hosted under a different GitHub username, please replace `zhenghaoshe` with the correct username.
