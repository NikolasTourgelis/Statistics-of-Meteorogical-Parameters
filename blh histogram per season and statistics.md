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

# Extract the 'blh' variable (3D data: time, lat, lon)
blh_2 = ds2.variables['blh'][:]

# Extract time series for Athens and Thessaloniki
blh_athens = blh_2[:, athens_lat_idx, athens_lon_idx]
blh_thessaloniki = blh_2[:, thessaloniki_lat_idx, thessaloniki_lon_idx]

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
seasons_athens, blh_stats_athens = compute_seasonal_stats(dates, blh_athens)
seasons_thessaloniki, blh_stats_thessaloniki = compute_seasonal_stats(dates, blh_thessaloniki)

# Display the statistics
print("Seasonal Statistics for Athens (Athina):")
for season, stats in blh_stats_athens.items():
    print(f"{season}: {stats}")

print("\nSeasonal Statistics for Thessaloniki:")
for season, stats in blh_stats_thessaloniki.items():
    print(f"{season}: {stats}")

# Plot the seasonal average blh data for both Athens and Thessaloniki
fig, ax = plt.subplots(figsize=(10, 6))

# Width of the bars
bar_width = 0.3
# Set positions for Athens and Thessaloniki bars
positions_athens = np.arange(len(seasons_athens)) - bar_width / 2
positions_thessaloniki = np.arange(len(seasons_thessaloniki)) + bar_width / 2

# Get seasonal averages for plotting
blh_seasonal_athens = [blh_stats_athens[season]['average'] for season in seasons_athens]
blh_seasonal_thessaloniki = [blh_stats_thessaloniki[season]['average'] for season in seasons_thessaloniki]

# Histogram for Athens
ax.bar(positions_athens, blh_seasonal_athens, width=bar_width, color='blue', alpha=0.7, label='Athens (Athina)')

# Histogram for Thessaloniki
ax.bar(positions_thessaloniki, blh_seasonal_thessaloniki, width=bar_width, color='red', alpha=0.7, label='Thessaloniki')

# Set x-ticks to the middle of each grouped bar position
ax.set_xticks(np.arange(len(seasons_athens)))
ax.set_xticklabels(seasons_athens)

# Add labels and title
ax.set_xlabel('Season')
ax.set_ylabel('Average BLH (Boundary Layer Height)')
ax.set_title('Seasonal Avg BLH (Boundary Layer Height) for Athens and Thessaloniki')
ax.legend()

# Display the plot
plt.tight_layout()
plt.show()

```

    Seasonal Statistics for Athens (Athina):
    Winter: {'average': 469.4931681875674, 'median': 400.34752148747884, 'max': 1475.156151990751, 'std_dev': 361.7825836255507, 'variance': 130886.63781477859, 'range': 1454.7655378483155}
    Spring: {'average': 458.05512975927866, 'median': 318.68698701364747, 'max': 1824.5319739154838, 'std_dev': 391.4665031445078, 'variance': 153246.02308418893, 'range': 1808.99248063288}
    Summer: {'average': 504.09967987691596, 'median': 320.3515873086876, 'max': 2192.646439161527, 'std_dev': 440.34731099014675, 'variance': 193905.75429625303, 'range': 2163.6950235017416}
    Autumn: {'average': 450.6009420429428, 'median': 332.81230951727457, 'max': 1895.3963864757684, 'std_dev': 399.8510591134445, 'variance': 159880.86947414326, 'range': 1874.7204122827548}
    
    Seasonal Statistics for Thessaloniki:
    Winter: {'average': 317.0723247274584, 'median': 205.39904693405128, 'max': 1800.6568496837635, 'std_dev': 327.1968543679668, 'variance': 107057.78150829246, 'range': 1786.639276670911}
    Spring: {'average': 402.16410702684897, 'median': 240.78369320576394, 'max': 2339.416625175647, 'std_dev': 386.62069778667313, 'variance': 149475.56395705402, 'range': 2315.9821704937094}
    Summer: {'average': 433.9618498367179, 'median': 287.2973814500315, 'max': 2189.887958672603, 'std_dev': 403.7277045980718, 'variance': 162996.0594600279, 'range': 2169.2119844795893}
    Autumn: {'average': 362.0460355778951, 'median': 200.07232598992255, 'max': 1871.1407821766106, 'std_dev': 371.5693049045257, 'variance': 138063.7483472324, 'range': 1855.6012888940068}
    


    
![png](output_0_1.png)
    



```python

```
