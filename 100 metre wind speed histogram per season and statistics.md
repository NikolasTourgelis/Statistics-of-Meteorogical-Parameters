```python
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from collections import defaultdict

# File paths for the two NetCDF files
file1 = 'athensnew.nc'
file2 = 'thessalonikinew.nc'

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

# Extract the 'u100' and 'v100' variables (3D data: time, lat, lon)
u100_2 = ds2.variables['u100'][:]
v100_2 = ds2.variables['v100'][:]

# Calculate wind speed from u100 and v100
wind_speed_2 = np.sqrt(u100_2**2 + v100_2**2)

# Extract time series for Athens and Thessaloniki
wind_speed_athens = wind_speed_2[:, athens_lat_idx, athens_lon_idx]
wind_speed_thessaloniki = wind_speed_2[:, thessaloniki_lat_idx, thessaloniki_lon_idx]

# Extract the 'valid_time' variable and convert it to datetime objects
valid_time = ds2.variables['valid_time'][:]
time_units = ds2.variables['valid_time'].units

# Convert valid_time to datetime using num2date from netCDF4
dates = nc.num2date(valid_time, units=time_units)

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
ax.set_ylabel('Average Wind Speed (100m)')
ax.set_title('Seasonal Avg Wind Speed (100m) for Athens and Thessaloniki')
ax.legend()

# Display the plot
plt.tight_layout()
plt.show()

```

    Seasonal Statistics for Athens (Athina):
    Winter: {'average': 3.2538717, 'median': 2.739933, 'max': 10.741056, 'std_dev': 2.0104141, 'variance': 4.041765, 'range': 10.635562}
    Spring: {'average': 3.1612675, 'median': 2.8864918, 'max': 9.202795, 'std_dev': 1.7399994, 'variance': 3.0275981, 'range': 9.136222}
    Summer: {'average': 3.1920745, 'median': 2.9768262, 'max': 8.205944, 'std_dev': 1.6196352, 'variance': 2.6232183, 'range': 8.110135}
    Autumn: {'average': 2.9927733, 'median': 2.5120797, 'max': 11.927088, 'std_dev': 1.7375481, 'variance': 3.0190732, 'range': 11.685983}
    
    Seasonal Statistics for Thessaloniki:
    Winter: {'average': 3.8121686, 'median': 3.133554, 'max': 12.76431, 'std_dev': 2.6430855, 'variance': 6.9859004, 'range': 12.630195}
    Spring: {'average': 3.3307104, 'median': 3.041054, 'max': 12.522216, 'std_dev': 2.073891, 'variance': 4.301024, 'range': 12.380643}
    Summer: {'average': 3.1382573, 'median': 2.81668, 'max': 9.932279, 'std_dev': 1.8205928, 'variance': 3.314558, 'range': 9.810999}
    Autumn: {'average': 3.358877, 'median': 2.8729267, 'max': 13.389228, 'std_dev': 2.2191715, 'variance': 4.9247227, 'range': 13.199022}
    


    
![png](output_0_1.png)
    



```python

```
