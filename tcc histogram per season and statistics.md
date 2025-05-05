```python
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from collections import defaultdict

# File paths for the two NetCDF files
file1 = 'athinacloudcover.nc'
file2 = 'thessalonikicloudcover.nc'

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

# Extract the 'tcc' variable (3D data: time, lat, lon)
tcc_2 = ds2.variables['tcc'][:]

# Extract time series for Athens and Thessaloniki
tcc_athens = tcc_2[:, athens_lat_idx, athens_lon_idx] * 100  # Convert to percentage
tcc_thessaloniki = tcc_2[:, thessaloniki_lat_idx, thessaloniki_lon_idx] * 100  # Convert to percentage

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
seasons_athens, tcc_stats_athens = compute_seasonal_stats(dates, tcc_athens)
seasons_thessaloniki, tcc_stats_thessaloniki = compute_seasonal_stats(dates, tcc_thessaloniki)

# Display the statistics
print("Seasonal Statistics for Athens (Athina):")
for season, stats in tcc_stats_athens.items():
    print(f"{season}: {stats}")

print("\nSeasonal Statistics for Thessaloniki:")
for season, stats in tcc_stats_thessaloniki.items():
    print(f"{season}: {stats}")

# Plot the seasonal average tcc data for both Athens and Thessaloniki
fig, ax = plt.subplots(figsize=(10, 6))

# Width of the bars
bar_width = 0.3
# Set positions for Athens and Thessaloniki bars
positions_athens = np.arange(len(seasons_athens)) - bar_width / 2
positions_thessaloniki = np.arange(len(seasons_thessaloniki)) + bar_width / 2

# Get seasonal averages for plotting
tcc_seasonal_athens = [tcc_stats_athens[season]['average'] for season in seasons_athens]
tcc_seasonal_thessaloniki = [tcc_stats_thessaloniki[season]['average'] for season in seasons_thessaloniki]

# Histogram for Athens
ax.bar(positions_athens, tcc_seasonal_athens, width=bar_width, color='blue', alpha=0.7, label='Athens (Athina)')

# Histogram for Thessaloniki
ax.bar(positions_thessaloniki, tcc_seasonal_thessaloniki, width=bar_width, color='red', alpha=0.7, label='Thessaloniki')

# Set x-ticks to the middle of each grouped bar position
ax.set_xticks(np.arange(len(seasons_athens)))
ax.set_xticklabels(seasons_athens)

# Add labels and title
ax.set_xlabel('Season')
ax.set_ylabel('Average TCC (Total Cloud Cover) (%)')
ax.set_title('Seasonal Avg TCC (Total Cloud Cover) for Athens and Thessaloniki')
ax.legend()

# Display the plot
plt.tight_layout()
plt.show()

```

    Seasonal Statistics for Athens (Athina):
    Winter: {'average': 45.5896, 'median': 42.100525, 'max': 100.0, 'std_dev': 37.387672, 'variance': 1397.8381, 'range': 100.0}
    Spring: {'average': 50.25674, 'median': 48.2666, 'max': 100.0, 'std_dev': 35.902622, 'variance': 1288.9982, 'range': 100.0}
    Summer: {'average': 19.811447, 'median': 5.0735474, 'max': 99.749756, 'std_dev': 27.898855, 'variance': 778.3461, 'range': 99.749756}
    Autumn: {'average': 39.22438, 'median': 28.590393, 'max': 100.0, 'std_dev': 36.361988, 'variance': 1322.1941, 'range': 100.0}
    
    Seasonal Statistics for Thessaloniki:
    Winter: {'average': 47.41089, 'median': 47.90039, 'max': 100.0, 'std_dev': 40.0306, 'variance': 1602.449, 'range': 100.0}
    Spring: {'average': 57.856148, 'median': 63.81378, 'max': 100.0, 'std_dev': 36.024494, 'variance': 1297.7643, 'range': 100.0}
    Summer: {'average': 27.319098, 'median': 13.449097, 'max': 100.0, 'std_dev': 32.359337, 'variance': 1047.1266, 'range': 100.0}
    Autumn: {'average': 45.762547, 'median': 40.73639, 'max': 100.0, 'std_dev': 38.306583, 'variance': 1467.3943, 'range': 100.0}
    


    
![png](output_0_1.png)
    



```python

```
