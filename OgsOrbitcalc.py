# Activate the conda environment before running this script: conda activate orbitAnalyisis

# Import required libraries
import numpy as np
import matplotlib.pyplot as plt
from skyfield.api import load, EarthSatellite, wgs84  # For orbital calculations and satellite data
from datetime import datetime, timedelta
from skyfield.api import utc

# Define the location of the Optical Ground Station (OGS) in the Netherlands using WGS84 coordinates
ogs = wgs84.latlon(52.21475, 4.42914)

# Load the timescale from Skyfield (it is used to generate time objects for calculations)
ts = load.timescale()

# TLE data for Sentinel-2B (real TLE, updated manually as needed)
#tle_name = "Sentinel-2B"
#line1 = "1 42063U 17036A   24049.50767361  .00000200  00000-0  10231-4 0  9998"
#line2 = "2 42063  97.4281  17.4623 0001447  92.3241 267.8097 14.30825964280312"


#line1 = "1 99999U 23001A   25150.00000000  .00000000  00000-0  00000-0 0  0001"
#line2 = "2 99999  97.7000  75.0000 0010000 000.0000 000.0000 15.30000000 00001"

tle_name ="Eagle-1"

line1 = "1 00001U 23001A   23001.00000000  .00000000  00000-0  00000-0 0  0001"
line2 = "2 00001  97.8000 295.0000 0001000  0.0000  0.0000 15.23000000 00001"






# Create the EarthSatellite object with the TLE data
satellite = EarthSatellite(line1, line2, tle_name, ts)

# Set dynamic elevation threshold (in degrees) for when additional calculations should occur.
high_elevation_threshold = 84  # Only passes where elevation exceeds this value are treated specially

# Define the simulation period: starting at the current UTC time for 20 days
start_time = ts.utc(2023, 1, 1)
end_time = start_time + timedelta(days=15)
num_days = (end_time.utc_datetime() - start_time.utc_datetime()).days  # Compute the number of days in simulation

# Time step settings:
# When the satellite is not in a pass (elevation <= 10°), use a normal time step in seconds.
# When the satellite is in a pass (elevation > 10°), use a finer time step in seconds.
normal_time_step = 180  # seconds (instead of 3 minutes)
pass_time_step = 1     # seconds

# Data storage lists to hold various calculated metrics for each pass
max_elevations = []              # Maximum elevation reached during each pass
pass_durations = []              # Duration (minutes) of each pass
pass_times = []                  # Time when maximum elevation was reached
elevation_threshold_durations = []  # Duration (minutes) when elevation goes above the high threshold
azimuth_changes = []             # Change in azimuth during the threshold period
azimuth_rates = []               # Rate of change of azimuth (°/s) during the threshold period
min_distances = []               # Minimum distance (km) between OGS and satellite during each pass

# Initialize pass tracking variables
max_elevation = 0
pass_start_time = None
pass_max_time = None
threshold_start_time = None
threshold_end_time = None
azimuth_start = None
azimuth_end = None
in_pass = False  # Flag to indicate if we are currently in a pass (satellite is above 10°)

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

            # If the pass reached above the high threshold, calculate additional metrics
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

            # Reset maximum elevation and pass flag for the next pass
            max_elevation = 0
            in_pass = False

        # Use a coarser time step when the satellite is not in a visible pass
        time_step = timedelta(seconds=normal_time_step)  

    # Increment time by the appropriate time step
    t = ts.utc(t.utc_datetime() + time_step)

# Compute additional statistics for reporting later
num_high_elevations = len([elev for elev in max_elevations if elev > high_elevation_threshold])
max_duration_above_threshold = max(elevation_threshold_durations) if elevation_threshold_durations else 0
min_duration_above_threshold = min(elevation_threshold_durations) if elevation_threshold_durations else 0
max_azimuth_change = max(azimuth_changes) if azimuth_changes else 0
max_azimuth_range = max(azimuth_changes) if azimuth_changes else 0
max_azimuth_rate = max(azimuth_rates) if azimuth_rates else 0

# Define a filename prefix for saving plots, includes satellite name, simulation days, and elevation threshold
filename_prefix = f"Eagle-1_{num_days}days_{high_elevation_threshold}deg"

# Generate and save multiple plots for various tracked data
plot_data = [
    ("max_elevations", max_elevations, "Maximum Elevation (°)", "Number of Passes"),
    ("pass_durations", pass_durations, "Pass Duration (minutes)", "Number of Passes"),
    ("duration_above_threshold", elevation_threshold_durations, f"Duration Above {high_elevation_threshold}° (minutes)", "Number of Passes"),
    ("azimuth_change", azimuth_changes, "Azimuth Change (°)", "Pass Index"),
    ("azimuth_rate", azimuth_rates, "Azimuth Rate (°/s)", "Pass Index")
]


plot_data = [
    ("max_elevations", max_elevations, "Maximum Elevation (°)", "Number of Passes"),
    ("pass_durations", [dur * 60 for dur in pass_durations], "Pass Duration (seconds)", "Pass Index"),
    ("azimuth_change", azimuth_changes, "Azimuth Change (°)", "Pass Index"),
    ("duration_above_threshold", [dur * 60 for dur in elevation_threshold_durations], f"Duration Above {high_elevation_threshold}° (seconds)", "Pass Index"),
    ("azimuth_rate", azimuth_rates, "Azimuth Rate (°/s)", "Pass Index")
]

# Loop through each type of data and create either a histogram or scatter plot based on the data type
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
    # Title includes the maximum value observed in the data during the simulation period
    plt.grid()
    plt.savefig(f"{filename_prefix}_{filename}.png")

# Plot 6: Scatter Plot of Maximum Elevation vs. Time
plt.figure(figsize=(10, 5))
# Color each point red if above the high elevation threshold, otherwise blue
colors = ['red' if elev > high_elevation_threshold else 'blue' for elev in max_elevations]
plt.scatter(pass_times, max_elevations, c=colors, marker='o')
plt.xlabel('Date')
plt.ylabel('Maximum Elevation (°)')
# Calculate the highest elevation among all passes for inclusion in the title
highest_elev = max(max_elevations) if max_elevations else 0
plt.title(f'Max Elevation Over {num_days} Days | Count > {high_elevation_threshold}°: {num_high_elevations} | Highest Elevation: {highest_elev:.2f}° | Total Passes: {len(max_elevations)}')
plt.xticks(rotation=45)
plt.grid()
plt.savefig(f"{filename_prefix}_max_elevation_vs_time.png")

# Plot 7: Scatter Plot of Duration Above Elevation Threshold vs. Pass Index (Seconds)
plt.figure(figsize=(8, 5))
plt.scatter(range(len(elevation_threshold_durations)), 
            [dur * 60 for dur in elevation_threshold_durations],  # Convert durations from minutes to seconds
            color='red', edgecolor='black')
plt.xlabel('Pass Index')
plt.ylabel(f'Duration Above {high_elevation_threshold}° (seconds)')
plt.title(f'Duration Above {high_elevation_threshold}° Over {num_days} Days')
plt.grid()
plt.savefig(f"{filename_prefix}_duration_above_threshold_scatter.png")

# New Plot: Scatter Plot of Minimum Distance (km) for Passes Above Elevation Threshold
plt.figure(figsize=(8, 5))
plt.scatter(range(len(min_distances)), min_distances, color='green', edgecolor='black')
plt.xlabel('Pass Index')
plt.ylabel('Minimum Distance (km)')
plt.title(f'Minimum Distance for Passes Above {high_elevation_threshold}° Over {num_days} Days')
plt.grid()
plt.savefig(f"{filename_prefix}_min_distance_scatter.png")



# Plot 8: CCDF of Maximum Elevations
plt.figure(figsize=(8, 5))

# Sort max_elevations in ascending order
sorted_max_elev = np.sort(max_elevations)

# Compute the CCDF (1 - CDF)
ccdf_values = 1.0 - np.arange(1, len(sorted_max_elev) + 1) / len(sorted_max_elev)

plt.plot(sorted_max_elev, ccdf_values, marker='o', linestyle='-', color='blue')
plt.xlabel('Maximum Elevation (°)')
plt.ylabel('CCDF')
plt.title(f'CCDF of Maximum Elevations Over {num_days} Days | Passes > {high_elevation_threshold}°: {num_high_elevations}')
plt.grid()

# Restrict x-axis to values > 80°
plt.xlim(left=80)
plt.ylim(bottom=0, top=1)

# Save the CCDF plot
plt.savefig(f"{filename_prefix}_ccdf_max_elevations.png")
print(f"CCDF plot saved as {filename_prefix}_ccdf_max_elevations.png")


# (A duplicate CCDF block, if needed, update similarly)
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
print(f"CCDF plot saved as {filename_prefix}_ccdf_max_elevations.png")


# New Plot: Cumulative Count of Passes vs. Maximum Elevation
plt.figure(figsize=(8, 5))
# For each sorted maximum elevation, compute the absolute number of passes (i.e. number of passes with max elevation >= that value)
abs_counts = np.arange(len(sorted_max_elev), 0, -1)
plt.plot(sorted_max_elev, abs_counts, marker='o', linestyle='-', color='green')
plt.xlabel('Maximum Elevation (°)')
plt.ylabel('Number of Passes')
plt.title(f'Cumulative Count of Maximum Elevations Over {num_days} Days | Passes > {high_elevation_threshold}°: {num_high_elevations}')
plt.grid()
plt.xlim(left=80)
plt.savefig(f"{filename_prefix}_cumulative_max_elevations.png")
print(f"Cumulative count plot saved as {filename_prefix}_cumulative_max_elevations.png")

print(f"All plots saved with prefix: {filename_prefix}")




# Show all generated plots. This will display all active figures.
plt.show()