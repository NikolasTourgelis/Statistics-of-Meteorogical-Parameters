```python
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from collections import defaultdict

# File paths for the two NetCDF files
file1 = 'isshfa.nc'
file2 = 'isshft.nc'

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

# Extract the 'ishf' variable (3D data: time, lat, lon)
ishf_2 = ds2.variables['ishf'][:]

# Extract time series for Athens and Thessaloniki
ishf_athens = ishf_2[:, athens_lat_idx, athens_lon_idx]
ishf_thessaloniki = ishf_2[:, thessaloniki_lat_idx, thessaloniki_lon_idx]

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
seasons_athens, ishf_stats_athens = compute_seasonal_stats(dates, ishf_athens)
seasons_thessaloniki, ishf_stats_thessaloniki = compute_seasonal_stats(dates, ishf_thessaloniki)

# Display the statistics
print("Seasonal Statistics for Athens (Athina):")
for season, stats in ishf_stats_athens.items():
    print(f"{season}: {stats}")

print("\nSeasonal Statistics for Thessaloniki:")
for season, stats in ishf_stats_thessaloniki.items():
    print(f"{season}: {stats}")

# Plot the seasonal average ishf data for both Athens and Thessaloniki
fig, ax = plt.subplots(figsize=(10, 6))

# Width of the bars
bar_width = 0.3
# Set positions for Athens and Thessaloniki bars
positions_athens = np.arange(len(seasons_athens)) - bar_width / 2
positions_thessaloniki = np.arange(len(seasons_thessaloniki)) + bar_width / 2

# Get seasonal averages for plotting
ishf_seasonal_athens = [ishf_stats_athens[season]['average'] for season in seasons_athens]
ishf_seasonal_thessaloniki = [ishf_stats_thessaloniki[season]['average'] for season in seasons_thessaloniki]

# Histogram for Athens
ax.bar(positions_athens, ishf_seasonal_athens, width=bar_width, color='blue', alpha=0.7, label='Athens (Athina)')

# Histogram for Thessaloniki
ax.bar(positions_thessaloniki, ishf_seasonal_thessaloniki, width=bar_width, color='red', alpha=0.7, label='Thessaloniki')

# Set x-ticks to the middle of each grouped bar position
ax.set_xticks(np.arange(len(seasons_athens)))
ax.set_xticklabels(seasons_athens)

# Set labels and title
ax.set_xlabel('Season')
ax.set_ylabel('Average ISHF (W/m^2)')
ax.set_title('Seasonal ISHF Comparison: Athens vs Thessaloniki')
ax.legend()

# Display the plot
plt.tight_layout()
plt.show()

```

    Seasonal Statistics for Athens (Athina):
    Winter: {'average': -20.159248, 'median': -14.2039795, 'max': 13.090576, 'std_dev': 22.63969, 'variance': 512.55554, 'range': 148.43726}
    Spring: {'average': -15.680948, 'median': -7.7269287, 'max': 18.807129, 'std_dev': 21.519957, 'variance': 463.10852, 'range': 108.77905}
    Summer: {'average': -23.608929, 'median': -7.765747, 'max': 23.003418, 'std_dev': 36.52215, 'variance': 1333.8672, 'range': 134.16211}
    Autumn: {'average': -15.198842, 'median': -8.057007, 'max': 14.317383, 'std_dev': 20.6201, 'variance': 425.18854, 'range': 125.014404}
    
    Seasonal Statistics for Thessaloniki:
    Winter: {'average': -9.107294, 'median': -5.178955, 'max': 64.08667, 'std_dev': 23.266993, 'variance': 541.35297, 'range': 201.2561}
    Spring: {'average': -20.68946, 'median': -3.7231445, 'max': 33.25244, 'std_dev': 34.695415, 'variance': 1203.7719, 'range': 198.05249}
    Summer: {'average': -39.404274, 'median': -8.165161, 'max': 42.63672, 'std_dev': 58.06206, 'variance': 3371.2031, 'range': 259.6316}
    Autumn: {'average': -19.248165, 'median': -5.359497, 'max': 41.888916, 'std_dev': 33.58408, 'variance': 1127.8905, 'range': 172.51172}
    


    
![png](output_0_1.png)
    



```python

```
