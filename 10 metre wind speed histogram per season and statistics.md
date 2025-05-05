```python
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from collections import defaultdict

# File paths for the two NetCDF files
file1 = 'Athina.nc'
file2 = 'Thessaloniki.nc'

# Coordinates for Athens and Thessaloniki
athens_coords = (37.9838, 23.7275)
thessaloniki_coords = (40.6401, 22.9444)

# Open the second NetCDF file
ds2 = nc.Dataset(file2, 'r')

# Extract latitude and longitude arrays from the file
latitudes = ds2.variables['latitude'][:]
longitudes = ds2.variables['longitude'][:]

# Function to find the nearest index for a given coordinate
def find_nearest_idx(array, value):
    return np.abs(array - value).argmin()

# Find the closest grid points for Athens and Thessaloniki
athens_lat_idx = find_nearest_idx(latitudes, athens_coords[0])
athens_lon_idx = find_nearest_idx(longitudes, athens_coords[1])

thessaloniki_lat_idx = find_nearest_idx(latitudes, thessaloniki_coords[0])
thessaloniki_lon_idx = find_nearest_idx(longitudes, thessaloniki_coords[1])

# Extract the 'u10' and 'v10' variables (3D data: time, lat, lon)
u10_2 = ds2.variables['u10'][:]
v10_2 = ds2.variables['v10'][:]

# Calculate wind speed from u10 and v10
wind_speed_2 = np.sqrt(u10_2**2 + v10_2**2)

# Extract time series for Athens and Thessaloniki
wind_speed_athens = wind_speed_2[:, athens_lat_idx, athens_lon_idx]
wind_speed_thessaloniki = wind_speed_2[:, thessaloniki_lat_idx, thessaloniki_lon_idx]

# Extract the time variable and convert it to datetime objects
time = ds2.variables['time'][:]
time_units = ds2.variables['time'].units

# Convert time to datetime using num2date from netCDF4
dates = nc.num2date(time, units=time_units)

# Convert cftime objects to standard datetime objects
dates = [datetime(dt.year, dt.month, dt.day) for dt in dates]

# Close the NetCDF file
ds2.close()

# Function to compute seasonal statistics across all years
def compute_seasonal_stats(dates, data):
    seasonal_data = defaultdict(list)
    for date, value in zip(dates, data):
        # Determine season based on month
        if date.month in [12, 1, 2]:  # Winter
            season = 'Winter'
        elif date.month in [3, 4, 5]:  # Spring
            season = 'Spring'
        elif date.month in [6, 7, 8]:  # Summer
            season = 'Summer'
        else:  # Autumn (9, 10, 11)
            season = 'Autumn'
        
        # Group data by season
        seasonal_data[season].append(value)
    
    # Compute statistics for each season
    seasonal_stats = {}
    for season, values in seasonal_data.items():
        values = np.array(values)
        seasonal_stats[season] = {
            'average': np.mean(values),
            'median': np.median(values),
            'max': np.max(values),
            'std_dev': np.std(values),
            'variance': np.var(values),
            'range': np.ptp(values)  # range = max - min
        }
    
    # Sort seasons for consistent plotting order
    sorted_seasons = ['Winter', 'Spring', 'Summer', 'Autumn']
    sorted_stats = {season: seasonal_stats[season] for season in sorted_seasons}
    
    return sorted_seasons, sorted_stats

# Calculate seasonal statistics for Athens and Thessaloniki
seasons_athens, wind_speed_stats_athens = compute_seasonal_stats(dates, wind_speed_athens)
seasons_thessaloniki, wind_speed_stats_thessaloniki = compute_seasonal_stats(dates, wind_speed_thessaloniki)

# Display the statistics
print("Seasonal Statistics for Athens (Athina):")
for season, stats in wind_speed_stats_athens.items():
    print(f"{season}: {stats}")

print("\nSeasonal Statistics for Thessaloniki:")
for season, stats in wind_speed_stats_thessaloniki.items():
    print(f"{season}: {stats}")

# Plot the seasonal average wind speed data for both Athens and Thessaloniki
fig, ax = plt.subplots(figsize=(10, 6))

# Width of the bars
bar_width = 0.3
# Set positions for Athens and Thessaloniki bars
positions_athens = np.arange(len(seasons_athens)) - bar_width / 2
positions_thessaloniki = np.arange(len(seasons_thessaloniki)) + bar_width / 2

# Get seasonal averages for plotting
wind_speed_seasonal_athens = [wind_speed_stats_athens[season]['average'] for season in seasons_athens]
wind_speed_seasonal_thessaloniki = [wind_speed_stats_thessaloniki[season]['average'] for season in seasons_thessaloniki]

# Histogram for Athens
ax.bar(positions_athens, wind_speed_seasonal_athens, width=bar_width, color='blue', alpha=0.7, label='Athens (Athina)')

# Histogram for Thessaloniki
ax.bar(positions_thessaloniki, wind_speed_seasonal_thessaloniki, width=bar_width, color='red', alpha=0.7, label='Thessaloniki')

# Set x-ticks to the middle of each grouped bar position
ax.set_xticks(np.arange(len(seasons_athens)))
ax.set_xticklabels(seasons_athens)

# Add labels and title
ax.set_xlabel('Season')
ax.set_ylabel('Average Wind Speed (10m)')
ax.set_title('Seasonal Avg Wind Speed (10m) for Athens and Thessaloniki')
ax.legend()

# Display the plot
plt.tight_layout()
plt.show()

```

    Seasonal Statistics for Athens (Athina):
    Winter: {'average': 2.5345341998519175, 'median': 2.261273014202582, 'max': 8.347225358071356, 'std_dev': 1.4754257670725595, 'variance': 2.1768811941416506, 'range': 8.251142488754223}
    Spring: {'average': 2.385009543129999, 'median': 2.199095198283375, 'max': 7.086106671880837, 'std_dev': 1.2272006605851165, 'variance': 1.5060214613405463, 'range': 6.944660787308342}
    Summer: {'average': 2.391195657752161, 'median': 2.2001271438214784, 'max': 6.002625468276837, 'std_dev': 1.10744829194503, 'variance': 1.2264417193319643, 'range': 5.739092138720206}
    Autumn: {'average': 2.3142470588557744, 'median': 2.0336660401873794, 'max': 9.113681964853837, 'std_dev': 1.3002367227828586, 'variance': 1.6906155352731083, 'range': 8.978848188698628}
    
    Seasonal Statistics for Thessaloniki:
    Winter: {'average': 2.6586975563088817, 'median': 2.1298526450745956, 'max': 9.38300896219671, 'std_dev': 1.748769070687798, 'variance': 3.058193262594265, 'range': 9.22991741141486}
    Spring: {'average': 2.396572832569703, 'median': 2.1685381694249304, 'max': 9.181646276450412, 'std_dev': 1.4719328936768261, 'variance': 2.1665864434878346, 'range': 9.157194993355711}
    Summer: {'average': 2.348799219342376, 'median': 2.096586495455881, 'max': 7.716325621738386, 'std_dev': 1.3725446424586834, 'variance': 1.8838787955420349, 'range': 7.55065521295408}
    Autumn: {'average': 2.4370733494559054, 'median': 2.062585859301177, 'max': 9.777373648794915, 'std_dev': 1.553366780733726, 'variance': 2.412948355487059, 'range': 9.418702376688172}
    


    
![png](output_0_1.png)
    



```python

```
