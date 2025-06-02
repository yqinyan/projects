import numpy as np
from scipy.interpolate import interp1d
import os


# files path definition
input_folder = "E:/2023_2025/M2/Projet/Abaqus_simulation/experimental_data"
output_folder = "E:/2023_2025/M2/Projet/Abaqus_simulation/interpolated_exp_data"
os.makedirs(output_folder, exist_ok=True)

def load_exp_data(i):
    data_path = os.path.join(input_folder, f"exp_deformation_{i:03}.csv")

    if not os.path.exists(data_path):
        print(f"File not found: {data_path}")
        return None
    
    data = np.loadtxt(data_path, delimiter=",", skiprows=1)

    return data


def spatial_interp(i, N):

    data = load_exp_data(i)
    if data is None:
        print(f"None file i={i}")
        return None
    x = data[:, 1]
    y = data[:, 2]
    M = len(x)
    
    dx = np.diff(x)
    dy = np.diff(y)
    ds = np.sqrt(dx**2 + dy**2) 
    t = np.zeros(M)
    t[1:] = np.cumsum(ds)         
    t /= t[-1]      # normalize t             
    
    t_new = np.linspace(0, 1, N)
    
    # interpolate for x and y
    f_x = interp1d(t, x, kind='cubic', fill_value="extrapolate")
    f_y = interp1d(t, y, kind='cubic', fill_value="extrapolate")
    x_new = f_x(t_new)
    y_new = f_y(t_new)
    
    return np.column_stack((x_new, y_new))


# Load original time for each file
time_original = np.loadtxt("E:/2023_2025/M2/Projet/Abaqus_simulation/time_abaqus.csv", skiprows=1)
# Uniform time nodes
N = 290  # Number of time nodes
time_uniform = np.linspace(time_original[0], time_original[-1], N)


# Load and organize data by time step
data_dict = {}
N_points = 8000

for i, t in enumerate(time_original):
    data = spatial_interp(i, N_points)
    if data is not None:  # Filter invalid data
        data_dict[t] = data

# Check for valid data
if not data_dict:
    raise ValueError("Error: data_dict is empty. Check experimental data files")

# Calculate maximum points (from valid data)
max_points = max(len(data) for data in data_dict.values())

# Initialize storage for interpolated data
interpolated_data = []

# Perform interpolation for each point across all time steps
for point_index in range(max_points):  # Using the correct range
    x_coords = []
    y_coords = []
    times_with_data = []

    # Collect x, y, and corresponding time for the current point
    for t, data in data_dict.items():
        if point_index < len(data):  # Only consider points that exist at this time step
            x_coords.append(data[point_index, 0])
            y_coords.append(data[point_index, 1])
            times_with_data.append(t)

    # If there are enough data points to interpolate, perform the interpolation
    if len(times_with_data) > 1:
        x_interp = interp1d(times_with_data, x_coords, kind='cubic', fill_value="extrapolate")(time_uniform)
        y_interp = interp1d(times_with_data, y_coords, kind='cubic', fill_value="extrapolate")(time_uniform)
    else:
        # If not enough data points, fill with NaN
        x_interp = np.full_like(time_uniform, np.nan)
        y_interp = np.full_like(time_uniform, np.nan)

    interpolated_data.append(np.column_stack((x_interp, y_interp)))

# Save results as CSV files
for i, t in enumerate(time_uniform):
    output_data = np.zeros((len(interpolated_data), 3))
    for j in range(len(interpolated_data)):
        output_data[j, 0] = j + 1  # Point number
        output_data[j, 1] = interpolated_data[j][i, 0]  # Interpolated x
        output_data[j, 2] = interpolated_data[j][i, 1]  # Interpolated y


    output_path = os.path.join(output_folder, f"image_{i:03}.csv")

    
    np.savetxt(output_path, output_data, delimiter=",", fmt="%.5f", header="num,x,y", comments="")
    print("Interpolation completed! Interpolated data saved to:", output_path)
