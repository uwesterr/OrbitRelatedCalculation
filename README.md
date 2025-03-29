# Orbit Analysis Tool

A Python-based tool for analyzing satellite orbits and ground station visibility windows.

## Overview

This tool simulates satellite passes over a ground station and generates detailed statistics and visualizations about the passes. It's particularly useful for:

- Analyzing maximum elevation angles of satellite passes
- Computing the duration of passes above specified elevation thresholds
- Calculating azimuth changes and rates during high-elevation passes
- Determining minimum distances between satellite and ground station
- Visualizing pass statistics with various plots

## Requirements

- Python 3.8+
- Required packages:
  - numpy
  - matplotlib
  - skyfield
  - pandas
  - pyarrow (for Feather file support)
  - plotly (optional, for interactive plots)
  - kaleido (optional, for saving Plotly plots)

Install in a conda environment with:

```bash
conda create -n orbitAnalyisis python=3.9
conda activate orbitAnalyisis
pip install numpy matplotlib skyfield pandas pyarrow plotly kaleido

```
## Usage

### Running a New Simulation

```bash
python OgsOrbitcalc.py
```

This will:

1. Run a simulation using the TLE data specified in the script  
2. Save the results to the default data files  
3. Generate and display analysis plots  

### Using Existing Data

To analyze previously saved data without running a new simulation:

```bash
python OgsOrbitcalc.py --use-existing --data-file=Bird12m1s.feather
```

### Command-Line Arguments

- `--use-existing`: Use existing data instead of running a new simulation
- `--data-file`: Path to the Feather file to read/write data (default: satellite_passes.feather)

## Output Files

The simulation produces two main data files:

### Main Pass Data File

**Default filename:** satellite_passes.feather

Contains one row per satellite pass with the following information:  

- `pass_time`: Timestamp when maximum elevation occurred  
- `max_elevation`: Maximum elevation angle (degrees)  
- `pass_duration_minutes`: Total pass duration (minutes)  
- `duration_above_threshold_minutes`: Duration above threshold (for high passes)  
- `azimuth_change`: Change in azimuth during high elevation period (degrees)  
- `azimuth_rate`: Rate of change of azimuth (degrees/second)  
- `min_distance_km`: Minimum distance between satellite and ground station (km)  
- `link_id`: Unique identifier for the pass  

### Detailed Link Data File

**Default filename:** satellite_passes_links.feather (derived from main filename)

Contains detailed timestep-by-timestep data during visible passes:  

- `time_tag`: Timestamp for each step  
- `azimuth`: Azimuth angle at that time (degrees)  
- `elevation`: Elevation angle at that time (degrees)  
- `distance_km`: Distance at that time (km)  
- `link_id`: Pass identifier (matches with the main file)  

### Visualization Files

The script generates several plot files saved with a prefix based on the satellite name, duration, and elevation threshold:

1. Distribution of maximum elevations
2. Pass durations vs. pass index
3. Azimuth change vs. pass index
4. Duration above threshold vs. pass index
5. Azimuth rate vs. pass index
6. Maximum elevation vs. time (with color coding for passes above threshold)
7. CCDF of maximum elevations
8. Cumulative count of passes vs. maximum elevation

Example filenames:  

- Bird-1_365days_84deg_max_elevations.png  
- Bird-1_365days_84deg_pass_durations.png   
- Bird-1_365days_84deg_max_elevation_vs_time.png  