# conda activate orbitAnalyisis

# Import necessary libraries
import numpy as np
import matplotlib.pyplot as plt
from skyfield.api import load, EarthSatellite, wgs84
from datetime import datetime, timedelta
from skyfield.api import utc

# Define the location of the Optical Ground Station (OGS) in Esslingen
ogs = wgs84.latlon(48.7408, 9.3101)

# Load the TLE data of the ISS
ts = load.timescale()
line1 = '1 25544U 98067A   20256.51661597  .00000700  00000-0  19456-4 0  9995'
line2 = '2 25544  51.6445  34.0150 0001745  94.2681 265.8660 15.49315599244307'
satellite = EarthSatellite(line1, line2, 'ISS (ZARYA)', ts)

# Set simulation period dynamically
start_time = ts.utc(datetime.utcnow().replace(tzinfo=utc))
end_time = start_time + timedelta(days=2)  # Set number of days dynamically
num_days = (end_time.utc_datetime() - start_time.utc_datetime()).days  # Compute dynamically

# Time step settings
normal_time_step = 3  # minutes when ISS is not visible
pass_time_step = 1     # seconds when ISS is above 10°

# Data storage
max_elevations = []  # Stores max elevation per pass
pass_durations = []   # Stores duration of each pass (in minutes)
pass_times = []       # Stores times of max elevations

# Initialize pass tracking
current_pass = []
max_elevation = 0
pass_start_time = None
pass_max_time = None
in_pass = False  # Flag to track if we are in an ISS pass

# Loop through time
t = start_time
while t < end_time:
    # Compute ISS position relative to OGS
    difference = satellite - ogs
    topocentric = difference.at(t)
    alt, az, distance = topocentric.altaz()

    elevation = alt.degrees  # Elevation in degrees

    if elevation > 10:  # ISS is visible
        if not in_pass:
            pass_start_time = t.utc_datetime()
            in_pass = True

        current_pass.append(elevation)

        if elevation > max_elevation:
            max_elevation = elevation
            pass_max_time = t.utc_datetime()  # Store time of max elevation

        time_step = timedelta(seconds=pass_time_step)  # 1s inside a pass
    else:
        if in_pass:  # If we were in a pass and now it ended
            pass_end_time = t.utc_datetime()
            duration = (pass_end_time - pass_start_time).total_seconds() / 60  # Convert to minutes

            # Store pass data
            max_elevations.append(max_elevation)
            pass_durations.append(duration)
            pass_times.append(pass_max_time)

            # Reset tracking variables
            current_pass = []
            max_elevation = 0
            in_pass = False

        time_step = timedelta(minutes=normal_time_step)  # 3min outside a pass

    # Move to the next time step
    t = ts.utc(t.utc_datetime() + time_step)

# Identify occasions where max elevation > 86°
high_elevation_times = [time for max_elev, time in zip(max_elevations, pass_times) if max_elev > 86]
high_elevation_values = [max_elev for max_elev in max_elevations if max_elev > 86]
num_high_elevations = len(high_elevation_values)  # Count of high elevation events

# Plot 1: Histogram of Max Elevations
plt.figure(figsize=(8, 5))
plt.hist(max_elevations, bins=10, color='blue', edgecolor='black')
plt.xlabel('Maximum Elevation (°)')
plt.ylabel('Number of Passes')
plt.title(f'Histogram of Maximum ISS Elevations Over {num_days} Days')
plt.grid()

# Plot 2: Histogram of Pass Durations
plt.figure(figsize=(8, 5))
plt.hist(pass_durations, bins=10, color='green', edgecolor='black')
plt.xlabel('Pass Duration (minutes)')
plt.ylabel('Number of Passes')
plt.title(f'Histogram of ISS Pass Durations Over {num_days} Days')
plt.grid()

# Plot 3: Max Elevation vs Time with high elevations marked
plt.figure(figsize=(10, 5))
plt.scatter(pass_times, max_elevations, color='blue', label='Max Elevations')
plt.scatter(high_elevation_times, high_elevation_values, color='red', label='Max Elev > 86°', zorder=3)
plt.xlabel('Date')
plt.ylabel('Maximum Elevation (°)')
plt.title(f'Maximum ISS Elevation Over {num_days} Days (Red: Max Elev > 86°) | Count: {num_high_elevations}')
plt.xticks(rotation=45)
plt.legend()
plt.grid()

plt.show()