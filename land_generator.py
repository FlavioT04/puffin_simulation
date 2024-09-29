
class LandGenerator:
    """
    A class to represent a landmass (e.g., island or land mass) in the simulation.
    Each location has a name and a position (x, y coordinates).
    """

    def __init__(self, name, type, detector_coords):
        """
        Initialize a new location with a name and coordinates.

        Args:
        - name (str): The name of the location (e.g., Island G).
        - type (str): The type of detector (emitter or receiver).
        - detector_coords (list): list consisting of tuples (x, y) that represent the coordinates.
        """
        self.name = name   
        self.type = type 
        self.coordinates = detector_coords

    def get_coordinates(self):
        """
        Get the coordinates of the location.

        Returns:
        - list: List containing coordinates (x, y).
        """
        return self.coordinates
    
    def get_type(self):
        """
        Get the type of detector

        Returns:
        - String; type of detector ('emitter' or 'receiver')
        """
        return self.type
    
    def __str__(self):
        """
        Return a string representation of the location, including its name and coordinates.

        Returns:
        - str: A string describing the location and its position.
        """
        return f"{self.name}"
    