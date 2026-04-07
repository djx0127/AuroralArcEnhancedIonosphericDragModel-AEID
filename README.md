# AuroralArcEnhancedIonosphericDragModel-AEID
## OverView
Through the analysis of 14 years of orbital data from 8 satellite formations, this work provides evidence from a dynamics perspective for the extreme auroral-arc charging phenomena observed by satellites including DMSP. And a mathematical Auroral-Arc Enhanced Ionospheric Drag (AEID) Model is established to describe the ionospheric drag enhancement.

## Repository Structure
```
AEID-Model/ 
│
├── README.md 
├── LICENSE
│
├── src/
│   ├── data_processing/
│   │   ├── main_TLEDataProcessing.m
│   │       ├── TLEDataHandle.m
│   │       ├── TLEelments_Mutation.m
│   │       └── LLAGenerate.m
│   │   └── LLA2MLAT.m
│   │
│   ├── data_analysis/
│   │   ├── CorrelationAnalysis_F107Focus.m
│   │   └── CorrelationAnalysis_DstFocus.m
│   │
│   └── model_optimization/
│       ├── main_ParamsOptimization_1stLayer.m
│       ├── main_ParamsOptimization_2ndLayer.m
│           └── ModelOptimizationComputation_1and2.m
│       └── main_ParamsOptimization_3rdLayer.m
│           └── ModelOptimizationComputation_3.m
│
└── data/
    ├── TLE_raw/                  # original TLE data from CelesTrak(https://celestrak.org/)
    └── TLE_IRI_dataset/          # Processed TLE-IRI dataset
```

## Requirements
- MATLAB (recommended R2021a or later)
- Aerospace Toolbox (recommended)
- IRI model support (external or integrated)

## Usage

### 1. Data Processing
Run:
```matlab
main_TLEDataProcessing.m
```

This script:
- Converts TLE coordinates
- Removes maneuver-contaminated segments
- Generates LLA trajectories
- Builds the TLE–IRI dataset


### 2. Data Analysis
Solar activity correlation:
```matlab
CorrelationAnalysis_F107Focus.m
```
Geomagnetic activity correlation:
```matlab
CorrelationAnalysis_DstFocus.m
```

### 3. Model Optimization
```matlab
ParamsOptimizationMain_1stLayer.m
ParamsOptimizationMain_2ndLayer.m
ParamsOptimizationMain_3rdLayer.m
```
Supporting subfunctions:
- Drag modeling
- Activation function construction
- Orbit propagation

## Dataset Description
### 1. Raw TLE Data
Source: Celestrak(https://celestrak.org/)
### 2. TLE–IRI Dataset
MATLAB .mat format
Includes:
- Orbit states in TLE-intervals and sub-intervals;
- Ionospheric parameters in sub-intervals from IRI-2020 model.

## External Dependencies
This repository uses the **AACGM-v2** library for geomagnetic coordinate conversion, which is required in: src/data_processing/LLA2MLAT.m
### AACGM-v2
- Description: Altitude Adjusted Corrected Geomagnetic Coordinates
- Usage: Conversion from geographic latitude/longitude (LLA) to magnetic latitude (MLAT)
Please install AACGM-v2 from the official repository:
https://github.com/aburrell/aacgmv2
### Notes
- This repository does NOT redistribute AACGM-v2 due to its license.
- Users must install and configure AACGM-v2 separately before running the code.

## Citation
If you use this repository, please cite:
Ding J., et al. 2026.
Auroral-arc Enhanced Ionospheric Drag Model via Satellite Orbit Decay Analysis across a Full Solar Cycle.
Space Weather, under review.

##  External Tool Citation
If you use this repository, please also cite the AACGM-v2 model:
Shepherd, S. G. (2014).
Altitude-adjusted corrected geomagnetic coordinates: Definition and functional approximations.
Journal of Geophysical Research: Space Physics, 119(9), 7501–7521.

## Contact
For questions or collaborations, please contact the author Jixin Ding(djx0127@buaa.edu.cn; djx0127@163.com).
