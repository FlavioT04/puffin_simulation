import numpy as np
import geopandas as gpd

class Puffin():
    '''
    Args:
        - position (tuple): current position
        - line (geopandas.GeoSeries): LineString representing a vector with given direction 
        - direction (float): direction of vector (line)
        - steps (int): number of steps left
        - crs (int): coordinate reference system
    '''
    def __init__(self, position, vector, direction, steps=100, crs=4269):

        self.position = position
        self.vector = vector
        self.direction = direction
        self.steps = steps
        self.crs = crs

    def set_position(self, new_pos):
        self.position = new_pos

    def set_direction(self, direction):
        self.direction = direction

    def set_vector(self, new_vector):
        self.vector = new_vector

    def set_steps(self, n):
        self.steps = n

    def get_position(self):
        return self.position

    def get_direction(self):
        return self.direction 

    def get_vector(self):
        return self.vector
    
    def get_steps(self):
        return self.steps
    
    def move(self):
        dx = np.cos(self.direction) * 0.00005         # change in x and y
        dy = np.sin(self.direction) * 0.00005
        self.vector = self.vector.translate(dx, dy)
        self.steps -= 1

        # update position
        new_x = self.position[0] + dx
        new_y = self.position[1] + dy
        self.position = (new_x, new_y)

        return self.vector
