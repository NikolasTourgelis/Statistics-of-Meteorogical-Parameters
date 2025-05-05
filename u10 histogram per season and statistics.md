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

# Extract the 'u10' variable (3D data: time, lat, lon)
u10_2 = ds2.variables['u10'][:]

# Extract time series for Athens and Thessaloniki
u10_athens = u10_2[:, athens_lat_idx, athens_lon_idx]
u10_thessaloniki = u10_2[:, thessaloniki_lat_idx, thessaloniki_lon_idx]

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
seasons_athens, u10_stats_athens = compute_seasonal_stats(dates, u10_athens)
seasons_thessaloniki, u10_stats_thessaloniki = compute_seasonal_stats(dates, u10_thessaloniki)

# Display the statistics
print("Seasonal Statistics for Athens (Athina):")
for season, stats in u10_stats_athens.items():
    print(f"{season}: {stats}")

print("\nSeasonal Statistics for Thessaloniki:")
for season, stats in u10_stats_thessaloniki.items():
    print(f"{season}: {stats}")

# Plot the seasonal average u10 data for both Athens and Thessaloniki
fig, ax = plt.subplots(figsize=(10, 6))

# Width of the bars
bar_width = 0.3
# Set positions for Athens and Thessaloniki bars
positions_athens = np.arange(len(seasons_athens)) - bar_width / 2
positions_thessaloniki = np.arange(len(seasons_thessaloniki)) + bar_width / 2

# Get seasonal averages for plotting
u10_seasonal_athens = [u10_stats_athens[season]['average'] for season in seasons_athens]
u10_seasonal_thessaloniki = [u10_stats_thessaloniki[season]['average'] for season in seasons_thessaloniki]

# Histogram for Athens
ax.bar(positions_athens, u10_seasonal_athens, width=bar_width, color='blue', alpha=0.7, label='Athens (Athina)')

# Histogram for Thessaloniki
ax.bar(positions_thessaloniki, u10_seasonal_thessaloniki, width=bar_width, color='red', alpha=0.7, label='Thessaloniki')

# Set x-ticks to the middle of each grouped bar position
ax.set_xticks(np.arange(len(seasons_athens)))
ax.set_xticklabels(seasons_athens)

# Add labels and title
ax.set_xlabel('Season')
ax.set_ylabel('Average U10 (10-meter U Wind Component)')
ax.set_title('Seasonal Avg U10 (10-meter U Wind Component) for Athens and Thessaloniki')
ax.legend()

# Display the plot
plt.tight_layout()
plt.show()

```

    Seasonal Statistics for Athens (Athina):
    Winter: {'average': 0.13791140203479013, 'median': -0.12336694307810273, 'max': 5.479610749592599, 'std_dev': 1.9422907902189128, 'variance': 3.772493513769209, 'range': 9.900743068496979}
    Spring: {'average': 0.619874211547194, 'median': 0.6383738355718237, 'max': 6.623794462527542, 'std_dev': 1.9971999232134734, 'variance': 3.9888075332839037, 'range': 10.912833004787426}
    Summer: {'average': 0.6962862223017948, 'median': 0.7739272111284373, 'max': 5.38651608795743, 'std_dev': 1.7237297254318793, 'variance': 2.971244166337462, 'range': 9.743488572491627}
    Autumn: {'average': 0.5396698542220026, 'median': 0.5433921199845909, 'max': 9.092690049811573, 'std_dev': 1.726738269843743, 'variance': 2.981625052542963, 'range': 12.95965752479309}
    
    Seasonal Statistics for Thessaloniki:
    Winter: {'average': -0.36856688961053713, 'median': -0.45674647460944956, 'max': 4.600243607930669, 'std_dev': 1.684791256602011, 'variance': 2.8385215783225832, 'range': 12.165207810974053}
    Spring: {'average': -0.11630074648586222, 'median': -0.18438168752817943, 'max': 4.394554727155668, 'std_dev': 1.639925025696855, 'variance': 2.6893540899068302, 'range': 8.928910283183903}
    Summer: {'average': 0.1608278262515701, 'median': 0.08326546467293017, 'max': 4.029724296423251, 'std_dev': 1.2157193233864314, 'variance': 1.4779734732551624, 'range': 7.1091612554099655}
    Autumn: {'average': -0.13722261944042047, 'median': -0.3812643165268805, 'max': 7.837799171688862, 'std_dev': 1.6742426643237873, 'variance': 2.803088499042014, 'range': 11.716088970382767}
    


    
![png](output_0_1.png)
    



```python

```
