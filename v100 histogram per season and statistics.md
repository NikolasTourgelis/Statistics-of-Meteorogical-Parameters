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

# Extract the 'v100' variable (3D data: time, lat, lon)
v100_2 = ds2.variables['v100'][:]

# Extract time series for Athens and Thessaloniki
v100_athens = v100_2[:, athens_lat_idx, athens_lon_idx]
v100_thessaloniki = v100_2[:, thessaloniki_lat_idx, thessaloniki_lon_idx]

# Extract the time variable and convert it to datetime objects
time = ds2.variables['valid_time'][:]
time_units = ds2.variables['valid_time'].units

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
seasons_athens, v100_stats_athens = compute_seasonal_stats(dates, v100_athens)
seasons_thessaloniki, v100_stats_thessaloniki = compute_seasonal_stats(dates, v100_thessaloniki)

# Display the statistics
print("Seasonal Statistics for Athens (Athina):")
for season, stats in v100_stats_athens.items():
    print(f"{season}: {stats}")

print("\nSeasonal Statistics for Thessaloniki:")
for season, stats in v100_stats_thessaloniki.items():
    print(f"{season}: {stats}")

# Plot the seasonal average v100 data for both Athens and Thessaloniki
fig, ax = plt.subplots(figsize=(10, 6))

# Width of the bars
bar_width = 0.3
# Set positions for Athens and Thessaloniki bars
positions_athens = np.arange(len(seasons_athens)) - bar_width / 2
positions_thessaloniki = np.arange(len(seasons_thessaloniki)) + bar_width / 2

# Get seasonal averages for plotting
v100_seasonal_athens = [v100_stats_athens[season]['average'] for season in seasons_athens]
v100_seasonal_thessaloniki = [v100_stats_thessaloniki[season]['average'] for season in seasons_thessaloniki]

# Histogram for Athens
ax.bar(positions_athens, v100_seasonal_athens, width=bar_width, color='blue', alpha=0.7, label='Athens (Athina)')

# Histogram for Thessaloniki
ax.bar(positions_thessaloniki, v100_seasonal_thessaloniki, width=bar_width, color='red', alpha=0.7, label='Thessaloniki')

# Set x-ticks to the middle of each grouped bar position
ax.set_xticks(np.arange(len(seasons_athens)))
ax.set_xticklabels(seasons_athens)

# Add labels and title
ax.set_xlabel('Season')
ax.set_ylabel('Average V100 (100-meter V Wind Component)')
ax.set_title('Seasonal Avg V100 (100-meter V Wind Component) for Athens and Thessaloniki')
ax.legend()

# Display the plot
plt.tight_layout()
plt.show()

```

    Seasonal Statistics for Athens (Athina):
    Winter: {'average': -0.6765646, 'median': -0.8463135, 'max': 8.776062, 'std_dev': 2.831072, 'variance': 8.014969, 'range': 18.986862}
    Spring: {'average': -1.063968, 'median': -1.1864471, 'max': 5.361862, 'std_dev': 1.9300231, 'variance': 3.724989, 'range': 12.5854645}
    Summer: {'average': -1.9151651, 'median': -1.8846817, 'max': 4.507431, 'std_dev': 1.749056, 'variance': 3.0591967, 'range': 12.171494}
    Autumn: {'average': -1.0873566, 'median': -1.1526642, 'max': 5.7480316, 'std_dev': 2.3417952, 'variance': 5.484005, 'range': 14.133926}
    
    Seasonal Statistics for Thessaloniki:
    Winter: {'average': -1.1590841, 'median': -1.4084549, 'max': 10.868622, 'std_dev': 3.8050225, 'variance': 14.478195, 'range': 22.311447}
    Spring: {'average': 0.036190696, 'median': 0.07406616, 'max': 6.577835, 'std_dev': 3.2027557, 'variance': 10.257645, 'range': 17.708176}
    Summer: {'average': 0.515406, 'median': 0.18269348, 'max': 8.366531, 'std_dev': 3.1982224, 'variance': 10.228627, 'range': 17.443848}
    Autumn: {'average': -0.25485578, 'median': -0.30343628, 'max': 7.7949066, 'std_dev': 3.3355718, 'variance': 11.1260395, 'range': 19.194733}
    


    
![png](output_0_1.png)
    



```python

```
