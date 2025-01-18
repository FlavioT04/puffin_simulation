import numpy as np
import geopandas as gpd
from shapely.geometry import LineString
from vector import Vector

class Puffin():
    """Puffin Class
    Description:
    ...

    Attributes:
        - position (geopandas.GeoSeries):
        - direction (angle):
        - vector (Numpy Array): 
        - bounces (int): 
    """
    def __init__(self, origin, direction, crs):
        self.position = origin
        self.direction = direction
        self.crs = crs
        self.vector = Vector(origin, direction, self.crs)
        self.bounces = 50

    
    def get_position(self):
        """ Returns the current position of the puffin """
        return self.position
    
    
    def get_direction(self):
        """ Returns the direction of the puffin """
        return self.direction
    

    def get_vector(self):
        v = self.vector
        return v.get_geoseries_vector()
    

    def get_bounces(self):
        return self.bounces
    

    def set_position(self, new_position):
        self.position = new_position


    def set_direction(self, new_direction):
        self.direction = new_direction


    def set_vector(self, new_vector):
        self.vector.set_geoseries_vector(new_vector)

    
    def decrease_bounces(self):
        self.bounces -= 1
    
    
    def move(self, object, direction):
        """
        Moves a puffin in the given direction

        Parameters:
            - object (geopandas.GeoSeries): Geometry representing the object
            - direction (float): angle of direction

        Returns:
            - translated_object (geopandas.GeoSeries): Translated geometry
        """
        dx = np.cos(direction) * 0.00006
        dy = np.sin(direction) * 0.00006
        new_position = gpd.GeoSeries.translate(object, dx, dy)
        return new_position
