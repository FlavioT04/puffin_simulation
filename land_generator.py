
class LandGenerator:
    """
    A class to represent a location (e.g., island or land mass) in the simulation.
    Each location has a name and a position (x, y coordinates).
    """

    def __init__(self, name, x, y):
        """
        Initialize a new location with a name and coordinates.

        Args:
        - name (str): The name of the location (e.g., Island G).
        - x (float): The x-coordinate of the location on the map.
        - y (float): The y-coordinate of the location on the map.
        """
        self.name = name    # Name of the location
        self.position = (x, y)  # Coordinates of the location

    def get_position(self):
        """
        Get the coordinates of the location.

        Returns:
        - tuple: The (x, y) coordinates of the location.
        """
        return self.position
    
    def __str__(self):
        """
        Return a string representation of the location, including its name and coordinates.

        Returns:
        - str: A string describing the location and its position.
        """
        return f"Location {self.name} at position {self.position}"
    