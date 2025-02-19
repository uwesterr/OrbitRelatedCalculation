# conda activate orbitAnalyisis

# Import necessary libraries
import numpy as np
import matplotlib.pyplot as plt
from skyfield.api import load, EarthSatellite, wgs84
from datetime import datetime, timedelta
from skyfield.api import utc

# Define the location of the Optical Ground Station (OGS) in the Netherlands
ogs = wgs84.latlon(52.21475, 4.42914)

# Load the timescale
ts = load.timescale()

# Using a real TLE for Sentinel-2B (Updated manually if needed)
tle_name = "Sentinel-2B"
line1 = "1 42063U 17036A   24049.50767361  .00000200  00000-0  10231-4 0  9998"
line2 = "2 42063  97.4281  17.4623 0001447  92.3241 267.8097 14.30825964280312"

# Create the EarthSatellite object
satellite = EarthSatellite(line1, line2, tle_name, ts)

# Set elevation threshold dynamically
high_elevation_threshold = 86  # Change this value as needed

# Set simulation period (e.g., 2 days)
start_time = ts.utc(datetime.utcnow().replace(tzinfo=utc))
end_time = start_time + timedelta(days=20)
num_days = (end_time.utc_datetime() - start_time.utc_datetime()).days  # Compute dynamically

# Time step settings
normal_time_step = 3  # minutes when satellite is not visible
pass_time_step = 1     # seconds when satellite is above 10°

# Data storage
max_elevations = []  
pass_durations = []   
pass_times = []       
elevation_threshold_durations = []  
azimuth_changes = []  
azimuth_rates = []  

# Initialize pass tracking
max_elevation = 0
pass_start_time = None
pass_max_time = None
threshold_start_time = None
threshold_end_time = None
azimuth_start = None
azimuth_end = None
in_pass = False  

# Loop through time
t = start_time
while t < end_time:
    difference = satellite - ogs
    topocentric = difference.at(t)
    alt, az, distance = topocentric.altaz()

    elevation = alt.degrees  
    azimuth = az.degrees  

    if elevation > 10:  
        if not in_pass:
            pass_start_time = t.utc_datetime()
            in_pass = True
            threshold_start_time = None  
            azimuth_start = None  

        if elevation > max_elevation:
            max_elevation = elevation
            pass_max_time = t.utc_datetime()  

        if elevation > high_elevation_threshold:  
            if threshold_start_time is None:
                threshold_start_time = t.utc_datetime()  
                azimuth_start = azimuth  
            threshold_end_time = t.utc_datetime()  
            azimuth_end = azimuth  

        time_step = timedelta(seconds=pass_time_step)  
    else:
        if in_pass:  
            pass_end_time = t.utc_datetime()
            duration = (pass_end_time - pass_start_time).total_seconds() / 60  

            max_elevations.append(max_elevation)
            pass_durations.append(duration)
            pass_times.append(pass_max_time)

            if max_elevation > high_elevation_threshold and threshold_start_time is not None:
                duration_above_threshold = (threshold_end_time - threshold_start_time).total_seconds() / 60
                elevation_threshold_durations.append(duration_above_threshold)

                azimuth_change = abs(azimuth_end - azimuth_start)  
                azimuth_changes.append(azimuth_change)

                azimuth_rate = azimuth_change / (duration_above_threshold * 60)  
                azimuth_rates.append(azimuth_rate)

            max_elevation = 0
            in_pass = False

        time_step = timedelta(minutes=normal_time_step)  

    t = ts.utc(t.utc_datetime() + time_step)

# Compute statistics
num_high_elevations = len([elev for elev in max_elevations if elev > high_elevation_threshold])
max_duration_above_threshold = max(elevation_threshold_durations) if elevation_threshold_durations else 0
min_duration_above_threshold = min(elevation_threshold_durations) if elevation_threshold_durations else 0
max_azimuth_change = max(azimuth_changes) if azimuth_changes else 0
max_azimuth_range = max(azimuth_changes) if azimuth_changes else 0
max_azimuth_rate = max(azimuth_rates) if azimuth_rates else 0

# Define filename prefix for saving plots
filename_prefix = f"Sentinel2_{num_days}days_{high_elevation_threshold}deg"

# Generate and save plots
plot_data = [
    ("max_elevations", max_elevations, "Maximum Elevation (°)", "Number of Passes"),
    ("pass_durations", pass_durations, "Pass Duration (minutes)", "Number of Passes"),
    ("duration_above_threshold", elevation_threshold_durations, f"Duration Above {high_elevation_threshold}° (minutes)", "Number of Passes"),
    ("azimuth_change", azimuth_changes, "Azimuth Change (°)", "Pass Index"),
    ("azimuth_rate", azimuth_rates, "Azimuth Rate (°/s)", "Pass Index")
]

for filename, data, xlabel, ylabel in plot_data:
    plt.figure(figsize=(8, 5))
    plt.hist(data, bins=10, color='blue', edgecolor='black') if "elevation" in filename else plt.scatter(range(len(data)), data, color='orange', edgecolor='black')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(f"{xlabel} Over {num_days} Days | Max: {max(data) if data else 0:.2f}")
    plt.grid()
    plt.savefig(f"{filename_prefix}_{filename}.png")

# Plot 6: Scatter Plot of Maximum Elevation vs. Time
plt.figure(figsize=(10, 5))
colors = ['red' if elev > high_elevation_threshold else 'blue' for elev in max_elevations]
plt.scatter(pass_times, max_elevations, c=colors, marker='o')
plt.xlabel('Date')
plt.ylabel('Maximum Elevation (°)')
plt.title(f'Max Elevation Over {num_days} Days | Count > {high_elevation_threshold}°: {num_high_elevations}')
plt.xticks(rotation=45)
plt.grid()
plt.savefig(f"{filename_prefix}_max_elevation_vs_time.png")

# Plot 7: Scatter Plot of Duration Above Elevation Threshold vs. Pass Index (Seconds)
plt.figure(figsize=(8, 5))
plt.scatter(range(len(elevation_threshold_durations)), 
            [dur * 60 for dur in elevation_threshold_durations],  # Convert minutes to seconds
            color='red', edgecolor='black')
plt.xlabel('Pass Index')
plt.ylabel(f'Duration Above {high_elevation_threshold}° (seconds)')
plt.title(f'Duration Above {high_elevation_threshold}° Over {num_days} Days')
plt.grid()
plt.savefig(f"{filename_prefix}_duration_above_threshold_scatter.png")

print(f"All plots saved with prefix: {filename_prefix}")

#show plots
plt.show()