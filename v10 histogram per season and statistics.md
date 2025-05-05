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

# Extract the 'v10' variable (3D data: time, lat, lon)
v10_2 = ds2.variables['v10'][:]

# Extract time series for Athens and Thessaloniki
v10_athens = v10_2[:, athens_lat_idx, athens_lon_idx]
v10_thessaloniki = v10_2[:, thessaloniki_lat_idx, thessaloniki_lon_idx]

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
seasons_athens, v10_stats_athens = compute_seasonal_stats(dates, v10_athens)
seasons_thessaloniki, v10_stats_thessaloniki = compute_seasonal_stats(dates, v10_thessaloniki)

# Display the statistics
print("Seasonal Statistics for Athens (Athina):")
for season, stats in v10_stats_athens.items():
    print(f"{season}: {stats}")

print("\nSeasonal Statistics for Thessaloniki:")
for season, stats in v10_stats_thessaloniki.items():
    print(f"{season}: {stats}")

# Plot the seasonal average v10 data for both Athens and Thessaloniki
fig, ax = plt.subplots(figsize=(10, 6))

# Width of the bars
bar_width = 0.3
# Set positions for Athens and Thessaloniki bars
positions_athens = np.arange(len(seasons_athens)) - bar_width / 2
positions_thessaloniki = np.arange(len(seasons_thessaloniki)) + bar_width / 2

# Get seasonal averages for plotting
v10_seasonal_athens = [v10_stats_athens[season]['average'] for season in seasons_athens]
v10_seasonal_thessaloniki = [v10_stats_thessaloniki[season]['average'] for season in seasons_thessaloniki]

# Histogram for Athens
ax.bar(positions_athens, v10_seasonal_athens, width=bar_width, color='blue', alpha=0.7, label='Athens (Athina)')

# Histogram for Thessaloniki
ax.bar(positions_thessaloniki, v10_seasonal_thessaloniki, width=bar_width, color='red', alpha=0.7, label='Thessaloniki')

# Set x-ticks to the middle of each grouped bar position
ax.set_xticks(np.arange(len(seasons_athens)))
ax.set_xticklabels(seasons_athens)

# Add labels and title
ax.set_xlabel('Season')
ax.set_ylabel('Average V10 (10-meter V wind Component)')
ax.set_title('Seasonal Avg V10 (10-meter V Wind Component) for Athens and Thessaloniki')
ax.legend()

# Display the plot
plt.tight_layout()
plt.show()

```

    Seasonal Statistics for Athens (Athina):
    Winter: {'average': -0.6051449372645893, 'median': -0.583708736418572, 'max': 5.831914841456605, 'std_dev': 2.1078499331506686, 'variance': 4.4430313406832775, 'range': 13.776117665679662}
    Spring: {'average': -0.8364185162837433, 'median': -0.8998451866820605, 'max': 3.6513990641728062, 'std_dev': 1.4565865838851653, 'variance': 2.121644476354256, 'range': 9.535996507201643}
    Summer: {'average': -1.3657236490910076, 'median': -1.3068118857152826, 'max': 3.039179595192394, 'std_dev': 1.273969637626239, 'variance': 1.622998637593531, 'range': 8.963883901314361}
    Autumn: {'average': -0.8382606058166534, 'median': -0.7848328586944108, 'max': 4.070161899409889, 'std_dev': 1.7523714039360054, 'variance': 3.070805537332647, 'range': 10.499740835057128}
    
    Seasonal Statistics for Thessaloniki:
    Winter: {'average': -0.8757581318946757, 'median': -1.0245893270378958, 'max': 7.245681765489408, 'std_dev': 2.5269646856451446, 'variance': 6.385550522497665, 'range': 15.933630977365933}
    Spring: {'average': 0.04546928109078843, 'median': 0.014945175927399523, 'max': 4.874658388513244, 'std_dev': 2.281490823175954, 'variance': 5.205200376236092, 'range': 13.257087672709755}
    Summer: {'average': 0.49468570111178123, 'median': 0.3458267964457147, 'max': 6.458879480691845, 'std_dev': 2.3774321362530975, 'variance': 5.6521835624889665, 'range': 13.300143569853851}
    Autumn: {'average': -0.17070137862600998, 'median': -0.23454310478427143, 'max': 5.263930883240674, 'std_dev': 2.345467407388634, 'variance': 5.501217359122362, 'range': 13.311349899247519}
    


    
![png](output_0_1.png)
    



```python

```
