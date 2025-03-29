import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from skyfield.api import load, EarthSatellite, wgs84
from datetime import datetime, timedelta
from skyfield.api import utc
import os

# Load the timescale
ts = load.timescale()

# Default TLE for Sentinel-2B
tle_name = "Sentinel-2B"
line1 = "1 42063U 17036A   24049.50767361  .00000200  00000-0  10231-4 0  9998"
line2 = "2 42063  97.4281  17.4623 0001447  92.3241 267.8097 14.30825964280312"
satellite = EarthSatellite(line1, line2, tle_name, ts)

# ---- Streamlit UI ----
st.title("Satellite Pass Analysis: Sentinel-2B")
st.sidebar.header("Simulation Settings")

# User input for OGS coordinates
ogs_lat = st.sidebar.number_input("OGS Latitude (°)", value=52.21475, format="%.5f")
ogs_lon = st.sidebar.number_input("OGS Longitude (°)", value=4.42914, format="%.5f")
ogs = wgs84.latlon(ogs_lat, ogs_lon)

# User input for elevation threshold
high_elevation_threshold = st.sidebar.slider("Elevation Threshold (°)", min_value=10, max_value=90, value=86)

# User input for number of simulation days
num_days = st.sidebar.number_input("Simulation Duration (days)", min_value=1, max_value=30, value=2)

# Set the simulation period
start_time = ts.utc(datetime.utcnow().replace(tzinfo=utc))
end_time = start_time + timedelta(days=num_days)

# Time step settings
normal_time_step = timedelta(minutes=3)  
pass_time_step = timedelta(seconds=1)  

# ---- Satellite Pass Analysis ----
max_elevations, pass_durations, pass_times = [], [], []
elevation_threshold_durations, azimuth_changes, azimuth_rates = [], [], []

# Initialize pass tracking
max_elevation, pass_start_time, pass_max_time = 0, None, None
threshold_start_time, threshold_end_time, azimuth_start, azimuth_end = None, None, None, None
in_pass = False

# Process passes
t = start_time
while t < end_time:
    difference = satellite - ogs
    topocentric = difference.at(t)
    alt, az, distance = topocentric.altaz()
    
    elevation, azimuth = alt.degrees, az.degrees

    if elevation > 10:
        if not in_pass:
            pass_start_time, in_pass = t.utc_datetime(), True
            threshold_start_time, azimuth_start = None, None
        
        if elevation > max_elevation:
            max_elevation, pass_max_time = elevation, t.utc_datetime()

        if elevation > high_elevation_threshold:
            if threshold_start_time is None:
                threshold_start_time, azimuth_start = t.utc_datetime(), azimuth
            threshold_end_time, azimuth_end = t.utc_datetime(), azimuth
        
        time_step = pass_time_step
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

            max_elevation, in_pass = 0, False

        time_step = normal_time_step

    t = ts.utc(t.utc_datetime() + time_step)

# ---- Plotting ----
plot_filenames = []
output_folder = "plots"
os.makedirs(output_folder, exist_ok=True)
filename_prefix = f"Sentinel2_{num_days}days_{high_elevation_threshold}deg"

def save_plot(data, title, xlabel, ylabel, filename, is_histogram=True):
    plt.figure(figsize=(8, 5))
    if is_histogram:
        plt.hist(data, bins=10, color='blue', edgecolor='black')
    else:
        plt.scatter(range(len(data)), data, color='orange', edgecolor='black')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid()
    filepath = os.path.join(output_folder, filename)
    plt.savefig(filepath)
    plot_filenames.append(filepath)

save_plot(max_elevations, f'Histogram of Max Elevations ({num_days} days)', 
          'Maximum Elevation (°)', 'Number of Passes', f"{filename_prefix}_max_elevations.png")

save_plot(pass_durations, f'Histogram of Pass Durations ({num_days} days)', 
          'Pass Duration (minutes)', 'Number of Passes', f"{filename_prefix}_pass_durations.png")

save_plot(elevation_threshold_durations, f'Duration Above {high_elevation_threshold}°', 
          f'Duration Above {high_elevation_threshold}° (minutes)', 'Number of Passes', 
          f"{filename_prefix}_duration_above_threshold.png")

save_plot(azimuth_changes, f'Azimuth Change ({num_days} days)', 
          'Azimuth Change (°)', 'Pass Index', f"{filename_prefix}_azimuth_change.png", is_histogram=False)

save_plot(azimuth_rates, f'Azimuth Rate ({num_days} days)', 
          'Azimuth Rate (°/s)', 'Pass Index', f"{filename_prefix}_azimuth_rate.png", is_histogram=False)

# Plot 6: Scatter Plot of Maximum Elevation vs. Time
plt.figure(figsize=(10, 5))
colors = ['red' if elev > high_elevation_threshold else 'blue' for elev in max_elevations]
plt.scatter(pass_times, max_elevations, c=colors, marker='o')
plt.xlabel('Date')
plt.ylabel('Maximum Elevation (°)')
plt.title(f'Max Elevation ({num_days} days)')
plt.xticks(rotation=45)
plt.grid()
filepath = os.path.join(output_folder, f"{filename_prefix}_max_elevation_vs_time.png")
plt.savefig(filepath)
plot_filenames.append(filepath)

# Plot 7: Scatter Plot of Duration Above Elevation Threshold vs. Pass Index (Seconds)
plt.figure(figsize=(8, 5))
plt.scatter(range(len(elevation_threshold_durations)), 
            [dur * 60 for dur in elevation_threshold_durations], 
            color='red', edgecolor='black')
plt.xlabel('Pass Index')
plt.ylabel(f'Duration Above {high_elevation_threshold}° (seconds)')
plt.title(f'Duration Above {high_elevation_threshold}° ({num_days} days)')
plt.grid()
filepath = os.path.join(output_folder, f"{filename_prefix}_duration_above_threshold_scatter.png")
plt.savefig(filepath)
plot_filenames.append(filepath)

# ---- Display Saved Plots ----
st.subheader("Generated Plots")
for filepath in plot_filenames:
    st.image(filepath, caption=os.path.basename(filepath))

st.success(f"All plots saved in `{output_folder}` folder.")