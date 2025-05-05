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

# Extract the 't2m' variable (3D data: time, lat, lon)
t2m_2 = ds2.variables['t2m'][:]

# Convert the temperature from Kelvin to Celsius
t2m_2_celsius = t2m_2 - 273.15

# Extract time series for Athens and Thessaloniki in Celsius
t2m_athens = t2m_2_celsius[:, athens_lat_idx, athens_lon_idx]
t2m_thessaloniki = t2m_2_celsius[:, thessaloniki_lat_idx, thessaloniki_lon_idx]

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
seasons_athens, t2m_stats_athens = compute_seasonal_stats(dates, t2m_athens)
seasons_thessaloniki, t2m_stats_thessaloniki = compute_seasonal_stats(dates, t2m_thessaloniki)

# Display the statistics
print("Seasonal Statistics for Athens (Athina):")
for season, stats in t2m_stats_athens.items():
    print(f"{season}: {stats}")

print("\nSeasonal Statistics for Thessaloniki:")
for season, stats in t2m_stats_thessaloniki.items():
    print(f"{season}: {stats}")

# Plot the seasonal average t2m data for both Athens and Thessaloniki
fig, ax = plt.subplots(figsize=(10, 6))

# Width of the bars
bar_width = 0.3
# Set positions for Athens and Thessaloniki bars
positions_athens = np.arange(len(seasons_athens)) - bar_width / 2
positions_thessaloniki = np.arange(len(seasons_thessaloniki)) + bar_width / 2

# Get seasonal averages for plotting
t2m_seasonal_athens = [t2m_stats_athens[season]['average'] for season in seasons_athens]
t2m_seasonal_thessaloniki = [t2m_stats_thessaloniki[season]['average'] for season in seasons_thessaloniki]

# Histogram for Athens
ax.bar(positions_athens, t2m_seasonal_athens, width=bar_width, color='blue', alpha=0.7, label='Athens (Athina)')

# Histogram for Thessaloniki
ax.bar(positions_thessaloniki, t2m_seasonal_thessaloniki, width=bar_width, color='red', alpha=0.7, label='Thessaloniki')

# Set x-ticks to the middle of each grouped bar position
ax.set_xticks(np.arange(len(seasons_athens)))
ax.set_xticklabels(seasons_athens)

# Add labels and title
ax.set_xlabel('Season')
ax.set_ylabel('Average T2M (°C)')  # Label in Celsius
ax.set_title('Seasonal Avg T2M (2-meter Temperature) for Athens and Thessaloniki')
ax.legend()

# Display the plot
plt.tight_layout()
plt.show()

```

    Seasonal Statistics for Athens (Athina):
    Winter: {'average': 10.539041265204423, 'median': 10.612715899396619, 'max': 18.810335130237377, 'std_dev': 3.459690955858736, 'variance': 11.969461510050733, 'range': 17.24249154230978}
    Spring: {'average': 14.030533969974543, 'median': 13.816890063677391, 'max': 22.402012641472197, 'std_dev': 3.3103152287395514, 'variance': 10.958186913624989, 'range': 16.17193144818259}
    Summer: {'average': 25.6786856060914, 'median': 25.82151856876979, 'max': 35.24310418242294, 'std_dev': 3.5440961363589, 'variance': 12.560617423754083, 'range': 18.33932304978981}
    Autumn: {'average': 19.502663686067358, 'median': 19.52153982028858, 'max': 28.723102347111876, 'std_dev': 3.512955949678658, 'variance': 12.340859504382683, 'range': 20.691552809629286}
    
    Seasonal Statistics for Thessaloniki:
    Winter: {'average': 9.27657996421535, 'median': 9.321194453677009, 'max': 18.639570943444085, 'std_dev': 4.085679410740588, 'variance': 16.692776247349556, 'range': 20.05916235963656}
    Spring: {'average': 14.968286256441147, 'median': 14.7138711767235, 'max': 24.72703272319626, 'std_dev': 3.832150256498137, 'variance': 14.685375588378736, 'range': 18.464112263215668}
    Summer: {'average': 26.703048908823057, 'median': 26.504106183561674, 'max': 37.909652637733586, 'std_dev': 3.7430475459180537, 'variance': 14.010404931003164, 'range': 20.02163176913251}
    Autumn: {'average': 19.67777894681469, 'median': 19.94845028727181, 'max': 31.697401644555555, 'std_dev': 5.034755658125483, 'variance': 25.34876453702656, 'range': 27.727600264370608}
    


    
![png](output_0_1.png)
    



```python

```
