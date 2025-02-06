
class Puffin():
    '''
    Args:
        - position (tuple): current position
        - line (geopandas.GeoSeries): LineString representing a vector with given direction 
        - direction (float): direction of vector (line)
        - steps (int): number of steps left
        - crs (int): coordinate reference system
    '''
    def __init__(self, position, line, direction, steps=100, crs=4269):

        self.position = position
        self.line = line
        self.direction = direction
        self.steps = steps
        self.crs = crs

    def set_position(self, new_pos):
        self.position = new_pos

    def set_direction(self, direction):
        self.direction = direction

    def set_line(self, line):
        self.line = line

    def set_steps(self, n):
        self.steps = n

    def get_position(self):
        return self.position

    def get_direction(self):
        return self.direction 

    def get_line(self):
        return self.line
    
    def get_steps(self):
        return self.steps
