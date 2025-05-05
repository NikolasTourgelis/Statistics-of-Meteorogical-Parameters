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

# Extract the 'u100' variable (3D data: time, lat, lon)
u100_2 = ds2.variables['u100'][:]

# Extract time series for Athens and Thessaloniki
u100_athens = u100_2[:, athens_lat_idx, athens_lon_idx]
u100_thessaloniki = u100_2[:, thessaloniki_lat_idx, thessaloniki_lon_idx]

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
seasons_athens, u100_stats_athens = compute_seasonal_stats(dates, u100_athens)
seasons_thessaloniki, u100_stats_thessaloniki = compute_seasonal_stats(dates, u100_thessaloniki)

# Display the statistics
print("Seasonal Statistics for Athens (Athina):")
for season, stats in u100_stats_athens.items():
    print(f"{season}: {stats}")

print("\nSeasonal Statistics for Thessaloniki:")
for season, stats in u100_stats_thessaloniki.items():
    print(f"{season}: {stats}")

# Plot the seasonal average u100 data for both Athens and Thessaloniki
fig, ax = plt.subplots(figsize=(10, 6))

# Width of the bars
bar_width = 0.3
# Set positions for Athens and Thessaloniki bars
positions_athens = np.arange(len(seasons_athens)) - bar_width / 2
positions_thessaloniki = np.arange(len(seasons_thessaloniki)) + bar_width / 2

# Get seasonal averages for plotting
u100_seasonal_athens = [u100_stats_athens[season]['average'] for season in seasons_athens]
u100_seasonal_thessaloniki = [u100_stats_thessaloniki[season]['average'] for season in seasons_thessaloniki]

# Histogram for Athens
ax.bar(positions_athens, u100_seasonal_athens, width=bar_width, color='blue', alpha=0.7, label='Athens (Athina)')

# Histogram for Thessaloniki
ax.bar(positions_thessaloniki, u100_seasonal_thessaloniki, width=bar_width, color='red', alpha=0.7, label='Thessaloniki')

# Set x-ticks to the middle of each grouped bar position
ax.set_xticks(np.arange(len(seasons_athens)))
ax.set_xticklabels(seasons_athens)

# Add labels and title
ax.set_xlabel('Season')
ax.set_ylabel('Average U100 (100-meter U Wind Component)')
ax.set_title('Seasonal Avg U100 (100-meter U Wind Component) for Athens and Thessaloniki')
ax.legend()

# Display the plot
plt.tight_layout()
plt.show()

```

    Seasonal Statistics for Athens (Athina):
    Winter: {'average': 0.13649987, 'median': -0.19908142, 'max': 7.157196, 'std_dev': 2.47752, 'variance': 6.138105, 'range': 12.247879}
    Spring: {'average': 0.787624, 'median': 0.48735046, 'max': 9.202164, 'std_dev': 2.7466052, 'variance': 7.5438404, 'range': 14.608398}
    Summer: {'average': 0.67141086, 'median': 0.61019135, 'max': 7.183838, 'std_dev': 2.3737547, 'variance': 5.6347117, 'range': 12.920776}
    Autumn: {'average': 0.6033001, 'median': 0.49510193, 'max': 11.898819, 'std_dev': 2.2238357, 'variance': 4.9454455, 'range': 16.513687}
    
    Seasonal Statistics for Thessaloniki:
    Winter: {'average': -0.30930904, 'median': -0.23648834, 'max': 5.9917145, 'std_dev': 2.366683, 'variance': 5.601188, 'range': 17.03476}
    Spring: {'average': -0.16011657, 'median': -0.27689362, 'max': 6.286209, 'std_dev': 2.260545, 'variance': 5.110064, 'range': 13.440598}
    Summer: {'average': 0.18493678, 'median': 0.034370422, 'max': 5.6156006, 'std_dev': 1.6231891, 'variance': 2.634743, 'range': 9.395798}
    Autumn: {'average': -0.12525043, 'median': -0.32073212, 'max': 10.875381, 'std_dev': 2.23609, 'variance': 5.0000987, 'range': 15.673874}
    


    
![png](output_0_1.png)
    



```python

```
