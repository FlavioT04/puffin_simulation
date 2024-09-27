
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
        - type (str): The type of detector.
        - detector_coords (list): list with x and y coordinates of the detectors.
        """
        self.name = name   
        self.type = type 
        self.coordinates = detector_coords

    def get_position(self):
        """
        Get the coordinates of the location.

        Returns:
        - tuple: The (x, y) coordinates of the location.
        """
        return self.coordinates
    
    def __str__(self):
        """
        Return a string representation of the location, including its name and coordinates.

        Returns:
        - str: A string describing the location and its position.
        """
        return f"{self.name} Type: {self.type}"
    