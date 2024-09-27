from land_generator import LandGenerator
import numpy as np
import matplotlib.pyplot as plt


def simulate_puffin_movement(max_steps=100, step_size=0.1):
    """
    Simulates the movement of a puffin from island G (origin) towards a random direction
    and visualizes the until it reaches a location west to the origin.

    Parameters:
    - max_steps: Maximum numberof steps the puffin can take.
    - step_size: Distance covered in each step.

    Returns:
    - positions: A list of positions (x, y) the puffin has moved through.
    """

    # Instance of the LandGenerator class
    gull_island = LandGenerator("Gull Island", 0, 0)

    # Starting position
    puffin_position = np.array(gull_island.get_position(), dtype=float)

    # Store puffin's positions for visualization
    positions = [puffin_position.copy()]

    # Choose a random angle in radians
    angle = np.random.uniform((np.pi)/2, (3 * np.pi)/2)

    # Simulate the puffin's random movement
    for _ in range(max_steps):
        # Calculate the new position
        puffin_position += np.array([step_size * np.cos(angle), step_size * np.sin(angle)])

        # Stop if the puffin is sufficiently west of island G (x-coord < -1)
        if puffin_position[0] < -1:
            break

        # Store the new position
        positions.append(puffin_position.copy())

    # Convert positions to NumPy array
    positions = np.array(positions)

    # Plot the simulation
    plot_puffin_movement(positions)

    return positions


def plot_puffin_movement(positions):
    """
    Plots the puffin's movement based on a list of positions.
    """

    # Plotting
    plt.figure(figsize=(10, 6))
    plt.plot(positions[:, 0], positions[:, 1], marker='o', label='Puffin Path')
    plt.scatter(0, 0, color='red', label='Gull Island (Origin)', zorder=5)

    # Add grid, labels, and legend
    plt.xlim(-2, 2) # x-axis limits
    plt.ylim(-4, 4) # y-axis limits
    plt.axvline(x=-1, color='gray', linestyle='--', label='Witless Bay (Destination)')
    plt.title('Puffin Simulation')
    plt.xlabel('x Position')
    plt.ylabel('y Position')
    plt.legend()
    plt.grid()
    plt.show()


if __name__ == '__main__':
    simulate_puffin_movement()
