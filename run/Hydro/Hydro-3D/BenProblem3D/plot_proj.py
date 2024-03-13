import yt
import sys
import re
"""
This script creates a projection plot of the density field along the z-axis.
Input: dataset number or None to get the last one from the OutputLog.
"""

# Function 
def get_data_dump_from_input():
    if len(sys.argv) < 2:
        return get_last_dataset()
    return sys.argv[1]

# Function to find the last dataset in the OutputLog
def get_last_dataset(output_log='OutputLog'):
    with open(output_log, 'r') as log_file:
        lines = log_file.readlines()
        last_line = lines[-1]  # Get the last line
        match = re.search(r'(DD\d{4})', last_line)  # Regex to find dataset pattern
        if match:
            return match.group(1)  # Return the last dataset number
        else:
            raise RuntimeError('No dataset found in OutputLog')

# Specify the dataset dump or get the last one from the OutputLog
dataset_dump = get_data_dump_from_input() 

# Load the dataset
ds = yt.load(f"./{dataset_dump}/noh3D_{dataset_dump[-4:]}")

# Set density field units
ds.quan(1.0, "msun/kpc**2")  # Ensure the units are properly interpreted by yt

# Create a projection plot of the density field along the z-axis
axis = 'x'
p = yt.ProjectionPlot(ds, axis, 'density', weight_field=None)
p.set_unit('density', 'msun/kpc**2')

# Timestamp the plot
p.annotate_timestamp(corner='lower_right', text_args={'color': 'white', 'fontsize': 18})

# Save the plot to a file with the dataset number
plot_filename = f"proj_density_DD{dataset_dump[-4:]}_{axis}.png"
p.save(plot_filename)
print(f"Plot saved to {plot_filename}")
