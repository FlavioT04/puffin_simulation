
class Puffin():

    def __init__(self, position, line, direction, crs):

        self.posiiton = position
        self.direction = direction
        self.crs = crs
        self.line = line

    def set_position(self, pos):
        self.position = pos

    def set_direction(self, direction):
        self.direction = direction

    def set_line(self, line):
        self.line = line

    def get_position(self):
        return self.posiiton

    def get_direction(self):
        return self.direction 

    def get_line(self):
        return self.line
