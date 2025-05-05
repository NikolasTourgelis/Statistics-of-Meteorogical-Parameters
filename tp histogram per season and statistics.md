```python
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from collections import defaultdict

# File paths for the two NetCDF files
file1 = 'athensnewtotal.nc'
file2 = 'thessalonikinewtotal.nc'

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

# Extract the 'tp' variable (3D data: time, lat, lon)
tp_2 = ds2.variables['tp'][:]

# Extract time series for Athens and Thessaloniki
tp_athens = tp_2[:, athens_lat_idx, athens_lon_idx]
tp_thessaloniki = tp_2[:, thessaloniki_lat_idx, thessaloniki_lon_idx]

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
seasons_athens, tp_stats_athens = compute_seasonal_stats(dates, tp_athens)
seasons_thessaloniki, tp_stats_thessaloniki = compute_seasonal_stats(dates, tp_thessaloniki)

# Display the statistics
print("Seasonal Statistics for Athens (Athina):")
for season, stats in tp_stats_athens.items():
    print(f"{season}: {stats}")

print("\nSeasonal Statistics for Thessaloniki:")
for season, stats in tp_stats_thessaloniki.items():
    print(f"{season}: {stats}")

# Plot the seasonal average tp data for both Athens and Thessaloniki
fig, ax = plt.subplots(figsize=(10, 6))

# Width of the bars
bar_width = 0.3
# Set positions for Athens and Thessaloniki bars
positions_athens = np.arange(len(seasons_athens)) - bar_width / 2
positions_thessaloniki = np.arange(len(seasons_thessaloniki)) + bar_width / 2

# Get seasonal averages for plotting
tp_seasonal_athens = [tp_stats_athens[season]['average'] for season in seasons_athens]
tp_seasonal_thessaloniki = [tp_stats_thessaloniki[season]['average'] for season in seasons_thessaloniki]

# Histogram for Athens
ax.bar(positions_athens, tp_seasonal_athens, width=bar_width, color='blue', alpha=0.7, label='Athens (Athina)')

# Histogram for Thessaloniki
ax.bar(positions_thessaloniki, tp_seasonal_thessaloniki, width=bar_width, color='red', alpha=0.7, label='Thessaloniki')

# Set x-ticks to the middle of each grouped bar position
ax.set_xticks(np.arange(len(seasons_athens)))
ax.set_xticklabels(seasons_athens)

# Add labels and title
ax.set_xlabel('Season')
ax.set_ylabel('Average TP (Total Precipitation)(m)')
ax.set_title('Seasonal Avg TP (Total Precipitation) for Athens and Thessaloniki')
ax.legend()

# Display the plot
plt.tight_layout()
plt.show()

```

    Seasonal Statistics for Athens (Athina):
    Winter: {'average': 8.636051e-05, 'median': 0.0, 'max': 0.0041337013, 'std_dev': 0.00032634332, 'variance': 1.0649996e-07, 'range': 0.0041337013}
    Spring: {'average': 5.7337078e-05, 'median': 0.0, 'max': 0.0016927719, 'std_dev': 0.0001940276, 'variance': 3.764671e-08, 'range': 0.0016927719}
    Summer: {'average': 3.5529551e-06, 'median': 0.0, 'max': 0.00036144257, 'std_dev': 2.3883196e-05, 'variance': 5.7040705e-10, 'range': 0.00036144257}
    Autumn: {'average': 7.138933e-05, 'median': 0.0, 'max': 0.0022668839, 'std_dev': 0.0002861797, 'variance': 8.1898825e-08, 'range': 0.0022668839}
    
    Seasonal Statistics for Thessaloniki:
    Winter: {'average': 8.197626e-05, 'median': 0.0, 'max': 0.0023140907, 'std_dev': 0.00028019954, 'variance': 7.851179e-08, 'range': 0.0023140907}
    Spring: {'average': 9.0698835e-05, 'median': 0.0, 'max': 0.0036125183, 'std_dev': 0.00030730144, 'variance': 9.443418e-08, 'range': 0.0036125183}
    Summer: {'average': 2.684671e-05, 'median': 0.0, 'max': 0.0018725395, 'std_dev': 0.00015453936, 'variance': 2.3882414e-08, 'range': 0.0018725395}
    Autumn: {'average': 5.3251184e-05, 'median': 0.0, 'max': 0.0025110245, 'std_dev': 0.00020264639, 'variance': 4.1065558e-08, 'range': 0.0025110245}
    


    
![png](output_0_1.png)
    



```python

```
