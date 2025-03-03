# Activate the conda environment before running this script: conda activate orbitAnalyisis

# author: Uwe Sterr ST2C
# date:Feb 2025
# parameters to set
# 1. TLE data for the satellite of interest
#  the variable is called line1 and line2
# 2. Time period for the analysis
#  the variable is called start_time and end_time
# 3. Location of the OGS
# the variable is called ogs
# 4. Elevation threshold for pass analysis
#  the variable is called high_elevation_threshold

# TODO :
# 1. Store the results in a file
# 2. read the data back from the file 
# 3. make a switch to switch between file reading or calculation
# 4. work again on the analysis of the data
# Import required libraries
import numpy as np
import matplotlib.pyplot as plt
from skyfield.api import load, EarthSatellite, wgs84  # For orbital calculations and satellite data
from datetime import datetime, timedelta
from skyfield.api import utc
import pandas as pd
import os
import argparse

# Add command-line argument parsing
parser = argparse.ArgumentParser(description='Satellite orbit analysis')
parser.add_argument('--use-existing', action='store_true',
                    help='Use existing data from feather file instead of running a new simulation')
parser.add_argument('--data-file', type=str, default='satellite_passes.feather',
                    help='Feather file to read/write pass data')
args = parser.parse_args()

# Function to run the simulation and return results as DataFrames
def run_simulation():
    print(f"Running new simulation for {tle_name}...")
    
    # Define the location of the Optical Ground Station (OGS) in the Netherlands using WGS84 coordinates
    ogs = wgs84.latlon(52.21475, 4.42914)

    # Load the timescale from Skyfield (it is used to generate time objects for calculations)
    ts = load.timescale()

    # Create the EarthSatellite object with the TLE data
    satellite = EarthSatellite(line1, line2, tle_name, ts)

    # Data storage lists to hold various calculated metrics for each pass
    max_elevations = []              # Maximum elevation reached during each pass
    pass_durations = []              # Duration (minutes) of each pass
    pass_times = []                  # Time when maximum elevation was reached
    elevation_threshold_durations = []  # Duration (minutes) when elevation goes above the high threshold
    azimuth_changes = []             # Change in azimuth during the threshold period
    azimuth_rates = []               # Rate of change of azimuth (°/s) during the threshold period
    min_distances = []               # Minimum distance (km) between OGS and satellite during each pass
    pass_indices = []                # Index of passes exceeding the threshold

    # Initialize pass tracking variables
    max_elevation = 0
    pass_start_time = None
    pass_max_time = None
    threshold_start_time = None
    threshold_end_time = None
    azimuth_start = None
    azimuth_end = None
    in_pass = False  # Flag to indicate if we are currently in a pass (satellite is above 10°)
    pass_idx = 0     # Counter for passes

    # Loop through the simulation period
    t = start_time
    while t < end_time:
        # Calculate the satellite's position relative to the OGS
        difference = satellite - ogs
        topocentric = difference.at(t)
        alt, az, distance = topocentric.altaz()

        # Extract the elevation (altitude) and azimuth angles in degrees
        elevation = alt.degrees  
        azimuth = az.degrees  
        current_distance = distance.km  # Get current distance (in kilometers) from OGS to satellite

        # If the satellite is visible (elevation > 20°), process the pass
        if elevation > 20:  
            if not in_pass:
                # If no pass is currently active, initialize a new pass
                pass_start_time = t.utc_datetime()
                in_pass = True
                threshold_start_time = None  
                azimuth_start = None  
                # Initialize the minimum distance for this pass to a very high value
                pass_min_distance = float('inf')

            # Update the minimum distance observed during the pass
            pass_min_distance = min(pass_min_distance, current_distance)

            # Update the maximum elevation reached during the pass
            if elevation > max_elevation:
                max_elevation = elevation
                pass_max_time = t.utc_datetime()  

            # If the elevation exceeds the high threshold, begin tracking the threshold period for additional stats
            if elevation > high_elevation_threshold:  
                if threshold_start_time is None:
                    # Mark the start time and initial azimuth when threshold is crossed
                    threshold_start_time = t.utc_datetime()  
                    azimuth_start = azimuth  
                # Always update the end time and ending azimuth for the threshold period
                threshold_end_time = t.utc_datetime()  
                azimuth_end = azimuth  

            # Use a fine time step during a pass
            time_step = timedelta(seconds=pass_time_step)  
        else:
            # When elevation falls below 10°, finish processing the current pass (if one is active)
            if in_pass:  
                pass_end_time = t.utc_datetime()
                # Calculate the overall duration of the pass in minutes
                duration = (pass_end_time - pass_start_time).total_seconds() / 60  
                max_elevations.append(max_elevation)
                pass_durations.append(duration)
                pass_times.append(pass_max_time)
                
                # For passes that exceed the threshold, record additional metrics
                if max_elevation > high_elevation_threshold and threshold_start_time is not None:
                    # Duration when elevation was above the high threshold (in minutes)
                    duration_above_threshold = (threshold_end_time - threshold_start_time).total_seconds() / 60
                    elevation_threshold_durations.append(duration_above_threshold)

                    # Calculate the change in azimuth during the period above the threshold
                    azimuth_change = abs(azimuth_end - azimuth_start)  
                    azimuth_changes.append(azimuth_change)

                    # Compute azimuth rate in degrees per second over the duration above threshold
                    azimuth_rate = azimuth_change / (duration_above_threshold * 60)  
                    azimuth_rates.append(azimuth_rate)

                    # Save the minimum distance observed during the entire pass
                    min_distances.append(pass_min_distance)
                    
                    # Save pass index for high elevation passes
                    pass_indices.append(pass_idx)

                # Reset maximum elevation and pass flag for the next pass
                max_elevation = 0
                in_pass = False
                pass_idx += 1  # Increment pass counter

            # Use a coarser time step when the satellite is not in a visible pass
            time_step = timedelta(seconds=normal_time_step)  

        # Increment time by the appropriate time step
        t = ts.utc(t.utc_datetime() + time_step)

    # Create DataFrames from the results
    # DataFrame for all passes
    all_passes_df = pd.DataFrame({
        'pass_time': pass_times,
        'max_elevation': max_elevations,
        'pass_duration_minutes': pass_durations
    })
    
    # DataFrame for high elevation passes
    high_passes_df = pd.DataFrame({
        'pass_index': pass_indices,
        'duration_above_threshold_minutes': elevation_threshold_durations,
        'azimuth_change': azimuth_changes,
        'azimuth_rate': azimuth_rates,
        'min_distance_km': min_distances
    })
    
    # Create a combined DataFrame with all data
    combined_df = pd.DataFrame({
        'pass_time': pass_times,
        'max_elevation': max_elevations,
        'pass_duration_minutes': pass_durations
    })
    
    # Add columns for high elevation data (will be NaN for passes that don't exceed threshold)
    combined_df['duration_above_threshold_minutes'] = np.nan
    combined_df['azimuth_change'] = np.nan
    combined_df['azimuth_rate'] = np.nan
    combined_df['min_distance_km'] = np.nan
    
    # Fill in high elevation data where applicable
    for idx, row in high_passes_df.iterrows():
        pass_idx = row['pass_index']
        combined_df.loc[pass_idx, 'duration_above_threshold_minutes'] = row['duration_above_threshold_minutes']
        combined_df.loc[pass_idx, 'azimuth_change'] = row['azimuth_change']
        combined_df.loc[pass_idx, 'azimuth_rate'] = row['azimuth_rate']
        combined_df.loc[pass_idx, 'min_distance_km'] = row['min_distance_km']
    
    return combined_df, all_passes_df, high_passes_df

# Function to create plots from the data in the DataFrames
def create_plots(all_passes_df, high_passes_df):
    # Compute additional statistics
    max_elevations = all_passes_df['max_elevation'].tolist()
    pass_durations = all_passes_df['pass_duration_minutes'].tolist()
    pass_times = all_passes_df['pass_time'].tolist()
    
    # Extract high pass data if available
    elevation_threshold_durations = []
    azimuth_changes = []
    azimuth_rates = []
    min_distances = []
    
    if not high_passes_df.empty:
        elevation_threshold_durations = high_passes_df['duration_above_threshold_minutes'].tolist()
        azimuth_changes = high_passes_df['azimuth_change'].tolist()
        azimuth_rates = high_passes_df['azimuth_rate'].tolist()
        min_distances = high_passes_df['min_distance_km'].tolist()
    
    # Calculate statistics for reports
    num_high_elevations = len(high_passes_df)
    max_duration_above_threshold = max(elevation_threshold_durations) if elevation_threshold_durations else 0
    min_duration_above_threshold = min(elevation_threshold_durations) if elevation_threshold_durations else 0
    max_azimuth_change = max(azimuth_changes) if azimuth_changes else 0
    max_azimuth_rate = max(azimuth_rates) if azimuth_rates else 0
    
    # Define a filename prefix for saving plots
    filename_prefix = f"{tle_name}_{num_days}days_{high_elevation_threshold}deg"
    
    # Create plots
    plot_data = [
        ("max_elevations", max_elevations, "Maximum Elevation (°)", "Number of Passes"),
        ("pass_durations", [dur * 60 for dur in pass_durations], "Pass Duration (seconds)", "Pass Index"),
        ("azimuth_change", azimuth_changes, "Azimuth Change (°)", "Pass Index"),
        ("duration_above_threshold", [dur * 60 for dur in elevation_threshold_durations], f"Duration Above {high_elevation_threshold}° (seconds)", "Pass Index"),
        ("azimuth_rate", azimuth_rates, "Azimuth Rate (°/s)", "Pass Index")
    ]

    # Loop through each type of data and create plots
    for filename, data, xlabel, ylabel in plot_data:
        plt.figure(figsize=(8, 5))
        # For the "azimuth_change" and "azimuth_rate" plots swap the labels
        if filename in ["azimuth_change", "azimuth_rate", "duration_above_threshold", "pass_durations"]:
            plot_xlabel = ylabel
            plot_ylabel = xlabel
            plt.title(f"{xlabel} Over {num_days} Days | Max: {max(data) if data else 0:.2f}")
            
        else:
            plot_xlabel = xlabel
            plot_ylabel = ylabel
            plt.title(f"{xlabel} Over {num_days} Days | Max: {max(data) if data else 0:.2f}")

        # For elevation data, use a histogram; otherwise, use a scatter plot
        if "elevation" in filename:
            plt.hist(data, bins=10, color='blue', edgecolor='black')
        else:
            plt.scatter(range(len(data)), data, color='orange', edgecolor='black')

        plt.xlabel(plot_xlabel)
        plt.ylabel(plot_ylabel)
        plt.grid()
        plt.savefig(f"{filename_prefix}_{filename}.png")

    # Plot 6: Maximum Elevation vs. Time
    plt.figure(figsize=(10, 5))
    colors = ['red' if elev > high_elevation_threshold else 'blue' for elev in max_elevations]
    plt.scatter(pass_times, max_elevations, c=colors, marker='o')
    plt.xlabel('Date')
    plt.ylabel('Maximum Elevation (°)')
    highest_elev = max(max_elevations) if max_elevations else 0
    plt.title(f'Max Elevation Over {num_days} Days | Count > {high_elevation_threshold}°: {num_high_elevations} | Highest Elevation: {highest_elev:.2f}° | Total Passes: {len(max_elevations)}')
    plt.xticks(rotation=45)
    plt.grid()
    plt.savefig(f"{filename_prefix}_max_elevation_vs_time.png")

    # More plots omitted for brevity - add these back as needed
    
    # Create CCDF plot
    plt.figure(figsize=(8, 5))
    sorted_max_elev = np.sort(max_elevations)
    ccdf_values = 1.0 - np.arange(1, len(sorted_max_elev) + 1) / len(sorted_max_elev)
    plt.plot(sorted_max_elev, ccdf_values, marker='o', linestyle='-', color='blue')
    plt.xlabel('Maximum Elevation (°)')
    plt.ylabel('CCDF')
    plt.title(f'CCDF of Maximum Elevations Over {num_days} Days | Passes > {high_elevation_threshold}°: {num_high_elevations}')
    plt.grid()
    plt.xlim(left=80)
    plt.ylim(bottom=0, top=1)
    plt.savefig(f"{filename_prefix}_ccdf_max_elevations.png")
    
    # Create Cumulative Count plot
    plt.figure(figsize=(8, 5))
    abs_counts = np.arange(len(sorted_max_elev), 0, -1)
    plt.plot(sorted_max_elev, abs_counts, marker='o', linestyle='-', color='green')
    plt.xlabel('Maximum Elevation (°)')
    plt.ylabel('Number of Passes')
    plt.title(f'Cumulative Count of Maximum Elevations Over {num_days} Days | Passes > {high_elevation_threshold}°: {num_high_elevations}')
    plt.grid()
    plt.xlim(left=80)
    plt.savefig(f"{filename_prefix}_cumulative_max_elevations.png")
    
    print(f"All plots saved with prefix: {filename_prefix}")
    return filename_prefix

# Main execution starts here
if __name__ == "__main__":
    # Satellite parameters
    tle_name ="Eagle-1"
    line1 = "1 00001U 23001A   23001.00000000  .00000000  00000-0  00000-0 0  0001"
    line2 = "2 00001  97.8000 295.0000 0001000  0.0000  0.0000 15.23000000 00001"

    # Simulation settings
    high_elevation_threshold = 84
    
    # Load the timescale from Skyfield (needed for time objects even when loading data)
    ts = load.timescale()
    
    # Define the simulation period
    start_time = ts.utc(2023, 1, 1)
    end_time = start_time + timedelta(days=365)
    num_days = (end_time.utc_datetime() - start_time.utc_datetime()).days
    
    # Time step settings
    normal_time_step = 10  # seconds 
    pass_time_step = 1     # seconds
    
    data_file = args.data_file
    
    # Either load existing data or run a new simulation
    if args.use_existing and os.path.exists(data_file):
        print(f"Loading existing data from {data_file}...")
        try:
            combined_df = pd.read_feather(data_file)
            print(f"Successfully loaded data with {len(combined_df)} passes")
            
            # Convert timestamp column back to datetime if needed
            if not pd.api.types.is_datetime64_dtype(combined_df['pass_time']):
                combined_df['pass_time'] = pd.to_datetime(combined_df['pass_time'])
            
            # Create separate DataFrames for all passes and high passes
            all_passes_df = combined_df[['pass_time', 'max_elevation', 'pass_duration_minutes']]
            high_passes_df = combined_df.dropna(subset=['duration_above_threshold_minutes']).copy()
            high_passes_df['pass_index'] = high_passes_df.index
            
        except Exception as e:
            print(f"Error loading data: {e}")
            print("Running new simulation instead...")
            combined_df, all_passes_df, high_passes_df = run_simulation()
            
            # Save results to a Feather file
            combined_df.reset_index().to_feather(data_file)
            print(f"Saved new simulation data to {data_file}")
    else:
        # Run a new simulation
        combined_df, all_passes_df, high_passes_df = run_simulation()
        
        # Save results to a Feather file
        combined_df.reset_index().to_feather(data_file)
        print(f"Saved simulation data to {data_file}")
    
    # Generate plots regardless of whether we loaded or simulated data
    filename_prefix = create_plots(all_passes_df, high_passes_df)
    
    # Display all plots
    plt.show()