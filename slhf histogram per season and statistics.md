```python
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
from collections import defaultdict

# File paths for the two NetCDF files
file1 = 'Athina.nc'
file2 = 'Thessaloniki.nc'

# Coordinates for Athens and Thessaloniki
athens_coords = (37.9838, 23.7275)
thessaloniki_coords = (40.6401, 22.9444)

# Open the second NetCDF file
ds2 = nc.Dataset(file2, 'r')

# Extract latitude and longitude arrays
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

# Extract the 'slhf' variable (3D data: time, lat, lon)
slhf_2 = ds2.variables['slhf'][:]
slhf_units = ds2.variables['slhf'].units.lower()  # Ensure case insensitivity

# Extract time variable and convert it to datetime objects
time = ds2.variables['time'][:]
time_units = ds2.variables['time'].units
dates = nc.num2date(time, units=time_units)
dates = [datetime(dt.year, dt.month, dt.day) for dt in dates]

# Convert slhf values to W/m² if necessary
if 'j' in slhf_units and '/' not in slhf_units:  # If in Joules (e.g., J/m² per time step)
    time_diff = np.diff(time).mean()  # Average time step in hours or seconds
    if 'hours' in time_units:
        slhf_2 /= (time_diff * 3600)  # Convert J/m² to W/m²
    elif 'seconds' in time_units:
        slhf_2 /= time_diff  # Already in seconds, just divide

# Extract SLHF time series for both cities
slhf_athens = slhf_2[:, athens_lat_idx, athens_lon_idx]
slhf_thessaloniki = slhf_2[:, thessaloniki_lat_idx, thessaloniki_lon_idx]

# Close the NetCDF file
ds2.close()

# Function to compute seasonal statistics
def compute_seasonal_stats(dates, data):
    seasonal_data = defaultdict(list)
    for date, value in zip(dates, data):
        season = ('Winter' if date.month in [12, 1, 2] else
                  'Spring' if date.month in [3, 4, 5] else
                  'Summer' if date.month in [6, 7, 8] else
                  'Autumn')
        seasonal_data[season].append(value)
    
    seasonal_stats = {season: {
        'average': np.mean(values),
        'median': np.median(values),
        'max': np.max(values),
        'std_dev': np.std(values),
        'variance': np.var(values),
        'range': np.ptp(values)
    } for season, values in seasonal_data.items()}
    
    sorted_seasons = ['Winter', 'Spring', 'Summer', 'Autumn']
    return sorted_seasons, {season: seasonal_stats[season] for season in sorted_seasons}

# Compute seasonal statistics for both cities
seasons_athens, slhf_stats_athens = compute_seasonal_stats(dates, slhf_athens)
seasons_thessaloniki, slhf_stats_thessaloniki = compute_seasonal_stats(dates, slhf_thessaloniki)

# Display results
print("Seasonal Statistics for Athens:")
for season, stats in slhf_stats_athens.items():
    print(f"{season}: {stats}")

print("\nSeasonal Statistics for Thessaloniki:")
for season, stats in slhf_stats_thessaloniki.items():
    print(f"{season}: {stats}")

# Plot the seasonal average SLHF for both cities
fig, ax = plt.subplots(figsize=(10, 6))
bar_width = 0.3
positions_athens = np.arange(len(seasons_athens)) - bar_width / 2
positions_thessaloniki = np.arange(len(seasons_thessaloniki)) + bar_width / 2

slhf_seasonal_athens = [slhf_stats_athens[season]['average'] for season in seasons_athens]
slhf_seasonal_thessaloniki = [slhf_stats_thessaloniki[season]['average'] for season in seasons_thessaloniki]

ax.bar(positions_athens, slhf_seasonal_athens, width=bar_width, color='blue', alpha=0.7, label='Athens')
ax.bar(positions_thessaloniki, slhf_seasonal_thessaloniki, width=bar_width, color='red', alpha=0.7, label='Thessaloniki')
ax.set_xticks(np.arange(len(seasons_athens)))
ax.set_xticklabels(seasons_athens)
ax.set_xlabel('Season')
ax.set_ylabel('Average SLHF (W/m²)')
ax.set_title('Seasonal Average SLHF for Athens and Thessaloniki')
ax.legend()
plt.tight_layout()
plt.show()

```

    Seasonal Statistics for Athens:
    Winter: {'average': -8.896172028792908, 'median': -7.574555235035665, 'max': -0.13731235351597812, 'std_dev': 6.014771288932568, 'variance': 36.17747365816754, 'range': 33.36138132131338}
    Spring: {'average': -6.06132295715591, 'median': -4.639472276407528, 'max': -0.1640169958194666, 'std_dev': 4.633720027024296, 'variance': 21.471361288846044, 'range': 23.330146594223656}
    Summer: {'average': -8.442726586760104, 'median': -7.223753342958112, 'max': -0.44562958738345143, 'std_dev': 5.484635003079358, 'variance': 30.081221117003313, 'range': 28.01074208159757}
    Autumn: {'average': -9.677092291041097, 'median': -8.040672627969172, 'max': -0.7393806527217492, 'std_dev': 7.415531712377041, 'variance': 54.99011057726956, 'range': 49.15596339644544}
    
    Seasonal Statistics for Thessaloniki:
    Winter: {'average': -7.008214510305616, 'median': -4.991488015862522, 'max': 0.3822143240244865, 'std_dev': 6.35006348800705, 'variance': 40.323306301720265, 'range': 32.62336211583534}
    Spring: {'average': -9.113990783476059, 'median': -5.628757889013785, 'max': -0.13245696400626372, 'std_dev': 8.261937394965706, 'variance': 68.25960951833272, 'range': 41.443177160248915}
    Summer: {'average': -11.302973215770976, 'median': -8.36598372512068, 'max': -1.17879340335177, 'std_dev': 8.029005871166525, 'variance': 64.46493527922654, 'range': 52.65669923295642}
    Autumn: {'average': -9.189512475796755, 'median': -7.229822579845269, 'max': -0.032921479056919146, 'std_dev': 8.133121627941563, 'variance': 66.1476674148908, 'range': 46.02909255218321}
    


    
![png](output_0_1.png)
    



```python

```
