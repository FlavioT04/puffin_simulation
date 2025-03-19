
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
