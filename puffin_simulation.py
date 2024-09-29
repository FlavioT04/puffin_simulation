from land_generator import LandGenerator
import numpy as np
import matplotlib.pyplot as plt
import math


def simulate_puffin_movement(source, destination, max_steps=10000000, step_size=0.001):
    """
    Simulates the movement of a puffin from island G (origin) towards a random direction
    and visualizes the until it reaches a location west to the origin.

    Args:
    - source (LandGenerator): The location where puffins are launched
    - destination (LandGenerator): The location where puffins end up
    - max_steps: Maximum number of steps the puffin can take.
    - step_size: Distance covered in each step.

    Returns:
    - positions: A list of positions (x, y) the puffin has moved through.
    """
    # Convert list of coords into a numpy array
    source_coordinates = np.array(source.get_coordinates())
    destination_coordinates = np.array(destination.get_coordinates())

    radius = calculate_average_distance(destination_coordinates)


    # Starting position
    random_index = np.random.randint(0, source_coordinates.shape[0])
    puffin_position = np.array(source_coordinates[random_index])

    # Store puffin's positions for visualization
    positions = [puffin_position.copy()]

    # Choose a random angle in radians
    angle = np.random.uniform((5 * np.pi)/6, (7 * np.pi)/6)

    # Simulate the puffin's random movement
    for _ in range(max_steps):
        # Calculate the new position
        puffin_position += np.array([step_size * np.cos(angle), step_size * np.sin(angle)])

        # Stop if the puffin reached a coordinate in the destination
        if detect(puffin_position, destination_coordinates, radius):
            break

        # Store the new position
        positions.append(puffin_position.copy())

    # Convert positions to NumPy array
    positions = np.array(positions)

    # Plot the simulation
    plot_puffin_movement(positions, source_coordinates, destination_coordinates)

    return positions


def plot_puffin_movement(positions, source_coords, destination_coords):
    """
    Plots the puffin's movement based on a list of positions.
    """
    # Split array into x and y coordinates
    x_source = source_coords[:, 0]
    y_source = source_coords[:, 1]

    x_destination = destination_coords[:, 0]
    y_destination = destination_coords[:, 1]

    # Plotting
    plt.scatter(positions[:, 0], positions[:, 1], color='orange')   # Puffin's path

    plt.scatter(x_source, y_source, color='blue')

    plt.scatter(x_destination, y_destination, color='green')

    # Add grid, labels, and legend
    plt.title('Puffin Simulation')
    plt.xlabel('x Position')
    plt.ylabel('y Position')
    plt.legend()
    plt.grid()
    plt.show()


def detect(current_position, receiver_coords, radius):
    x_coords = [x for x, y in receiver_coords]
    y_coords = [y for x, y in receiver_coords]
    detected = False

    for i in range(0, len(receiver_coords)):
        x_r, y_r = x_coords[i], y_coords[i]

        # Calculate the distance from the detection point to current position
        distance = np.sqrt((current_position[0] - x_r)**2 + (current_position[1] - y_r)**2)

        if distance <= radius:
            detected = True
            return detected

    return detected


def calculate_average_distance(coordinates):
    """
    Calculates the average distance between all pairs of points.

    args:
    - coordinates (np.ndarray): A 2D NumPy array of shape (n, 2), where n is the number of points.
                                Each row represents a point with its x and y coordinates.

    Returns:
    - float: the average distance between all pairs of points.
    """
    distances = []
    num_points = len(coordinates)
    
    for i in range(num_points):
        for j in range(i + 1, num_points):  # Only calculate each pair once
            distance = np.sqrt((coordinates[i][0] - coordinates[j][0])**2 + (coordinates[i][1] - coordinates[j][1])**2)
            distances.append(distance)
    
    # Calculate the average distance
    average_distance = np.mean(distances)
    return average_distance/50


