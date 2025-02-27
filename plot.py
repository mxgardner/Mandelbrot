import matplotlib.pyplot as plt
import numpy as np

# Example data for Experiment A and Experiment B
threads = [1, 2, 3, 4, 5, 10, 50]

# Replace these values with the actual execution times for your experiments
exec_time_A = [3.444, 3.284, 2.805, 2.494, 2.452, 2.116, 1.912]  # Experiment A times
exec_time_B = [8.395, 6.376, 5.805, 5.384, 5.198, 4.746, 4.710]  # Experiment B times

# Create the plot
plt.figure(figsize=(10, 6))

# Plot data for Experiment A
plt.plot(threads, exec_time_A, label='Experiment A (mandel -x -.5 -y .5 -s 1 -m 2000)', marker='o', color='blue')

# Plot data for Experiment B
plt.plot(threads, exec_time_B, label='Experiment B (mandel -x 0.2869325 -y 0.0142905 -s .000001 -W 1024 -H 1024 -m 1000)', marker='s', color='red')

# Labels and title
plt.xlabel('Number of Threads')
plt.ylabel('Execution Time (seconds)')
plt.title('Execution Time vs Number of Threads')

# Log scale if needed
# plt.yscale('log')

# Show legend
plt.legend()

# Show the plot
plt.grid(True)
plt.show()
