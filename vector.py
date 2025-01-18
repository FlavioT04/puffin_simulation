import numpy as np
import geopandas as gpd
from shapely.geometry import LineString


class Vector():
    """
    A class to represent a two-dimensional vector using a LineString shapely object as the geometry of the vector and a geopandas.GeoSeries
    object to store the geometry.

    Attributes:
        - origin (geopandas.GeoSeries): A GeoSeries object storing the geometry of the point (x, y) for the origin of the vector.
        - direction (float): An angle in radians representing the direction of the vector.
        - crs (int): Coordinate Reference System needed to create the vector with geopandas.
        - vector (geopandas.GeoSeries): A GeoSeries object representing the vector

    Methods:
        - create_vector():
            Creates vector using a LineString as its geometry and a GeoSeries object to store the geometry.

        - get_geoseries_vector():
            Returns self.vector.
        
        - set_geoseries_vector():
            Updates the self.vector attribute.
    """
    def __init__(self, origin, direction, crs):
        self.origin = origin
        self.direction = direction
        self.crs = crs
        self.vector = self.create_vector()


    def create_vector(self):
        """
        Creates a vector by using a LineString as its geometry and a GeoSeries object to store the geometry,
        """
        coords_1 = self.origin
        dx = np.cos(self.direction) * 0.00006
        dy = np.sin(self.direction) * 0.00006
        coords_2 = gpd.GeoSeries.translate(self.origin, dx, dy)
        x_1 = coords_1.x.iloc[0] 
        y_1 = coords_1.y.iloc[0] 
        x_2 = coords_2.x.iloc[0]
        y_2 = coords_2.y.iloc[0]
        return gpd.GeoSeries(LineString([(x_1, y_1), (x_2, y_2)]), crs=self.crs)
    

    def get_geoseries_vector(self):
        """ Returns vector """
        return self.vector
    

    def set_geoseries_vector(self, v):
        """ 
        Updates vector to a new vector v.

        Parameters:
            - v (geopandas.GeoSeries): new vector
        """
        self.vector = v
