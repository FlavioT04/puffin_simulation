import pandas as pd
import geopandas as gpd
import numpy as np
from puffin import Puffin
from shapely import Point
import matplotlib.pyplot as plt


class Simulation():
    """Simulation Class
    Description:
        ...

    Attributes:
        - puffins (list): List containing the geometry of the puffins (empty when simluation starts)
        - source (geopandas.GeoDataFrame): GeoDataFrame containing the geometry of the source landmass
        - destination (geopandas.GeoDataFrame): GeoDataFrame containing the geometry of the destination landmass
        - num_of_puffins (int): Number of puffins to be used in the simulation
        - crs (int): Coordinate Reference System to be used in the simulation

    Methods:
        - start(): 
            Starts simulation

        - get_random_coordinate(space):
            Returns a random coordinate within the geometry of a given space

        - spawn_puffin(origin):
            Spawns a puffin within the geometry of the origin

        - random_angle(min, max):
            Returns a random angle within the interval [min, max] in radians
    """
    def __init__(self, source, destination, num_of_puffins, crs):
        self.puffins = []
        self.source = source
        self.destination = destination
        self.num_of_puffins = num_of_puffins
        self.crs = crs


    def start(self):
        """
        Starts simulation 

        Returns:
            - puffins (list): List containing the geometry of the puffins
        """
        print('Initializing simulation...\n')
        puffins_reached_destination = 0
        for i in range(0, self.num_of_puffins):
            time = 0
            initial_position = self.spawn_puffin(self.source)    # Puffin's initial position
            direction = self.random_angle((5*np.pi)/6, (7*np.pi)/6)    # Angle within [150, 210)
            puffin = Puffin(initial_position, direction, self.crs)    # Initialize instance of Puffin class
            p_vector = puffin.get_vector()
            destination_reached = False
            bounce = False
            while not destination_reached and puffin.get_bounces() > 0:
                time += 1
                if time == 10000:
                    break
                p_vector = puffin.move(puffin.get_vector(), puffin.get_direction())
                puffin.set_vector(p_vector)    # Updates the vector of the puffin
                puffin_coordinates = puffin.get_position()
                if False: 
                    fig, ax = plt.subplots(figsize=(12, 8))
                    ax.set_xlim([-52.83935, -52.74986])
                    ax.set_ylim([47.24039, 47.29664]) 
                    self.destination.plot(ax=ax, color='black', linewidth=0.1)
                    #print('direction:', puffin.get_direction())
                    puffin.get_vector().plot(ax=ax, color='blue', marker='o', linewidth=1, label='Points')
                    plt.show()
                if not self.within_limits(puffin_coordinates):
                    print('puffin killed')
                    break
                elif p_vector.intersects(self.destination.geometry).iloc[0]:
                    intersection_point = p_vector.intersection(self.destination.geometry)
                    intersection_line = self.get_line(intersection_point, self.destination.geometry)    # Line intersecting with vector
                    puffin.set_position(intersection_point)      # Updates the position of the puffin
                    angle = self.compute_angle(puffin.get_vector(), intersection_line, puffin.get_direction())
                    print('intersection angle:', angle)
                    if angle > np.pi/4 or puffin.get_bounces() == 0:
                        puffins_reached_destination +=1
                        self.puffins.append(puffin.get_vector().iloc[0])
                        destination_reached = True
                        print(f'Puffin {puffins_reached_destination} reached destination')
                        print()
                    else:
                        bounce = True
                        print('Old direction:', puffin.get_direction())
                        puffin.decrease_bounces()
                        x = puffin.get_vector().get_coordinates().iloc[0].x
                        y = puffin.get_vector().get_coordinates().iloc[0].y
                        if puffin.get_direction() >= np.pi:
                            new_vector = puffin.get_vector().rotate(angle=-(angle), use_radians=True, origin=(x, y))
                            new_direction = puffin.get_direction() - angle
                        else:
                            new_vector = puffin.get_vector().rotate(angle=(angle), use_radians=True, origin=(x, y))
                            new_direction = puffin.get_direction() + angle
                        puffin.set_vector(new_vector)
                        puffin.set_direction(new_direction)
                        print('New direction:', puffin.get_direction())
                        # puffin.get_vector().plot(ax=ax, color='red', marker='o', markersize=10, label='Points')
                        print('Puffin bouncing')

        print('Number of puffins:', len(self.puffins))
        print()
        return self.puffins


    def get_random_coordinate(self, geometry):
        """
        Returns a GeoSeries Point with coordinates within the geometry of the given space

        Parameters:
            - geometry (geopandas.GeoSeries): Geometry used to get a coordinate from it

        Returns:
            - point (geopandas.GeoSeries): Point with coordinates (x, y) obtained from geometry
        """
        random_point = geometry.sample_points(size=1)
        return random_point


    def spawn_puffin(self, origin):
        """
        Spawns a puffin within the specified origin

        MIGHT DELETE FUNCTION TO USE GET_RANDOM_COORDINATE INSTEAD!!!!!

        Parameters:
            - origin (geopandas.GeoDataFrame): GeoDataFrame containing geometry where puffin is to be spawned within

        Returns:
            - new_puffin (shapely.Point): Point with coordinates (x, y) within the origin
        """
        new_puffin = self.get_random_coordinate(origin)
        return new_puffin

    
    def random_angle(self, min, max):
        """
        Returns an angle within the half-open interval [min, max) in radians

        Parameters:
            - min (float): lowerbound of the interval
            - max (float): upperbound of the interval

        Returns:
            - angle (float): angle in radians
        """
        angle = np.random.uniform(min, max)
        return angle
    

    def within_limits(self, coordinates):
        xmin = -52.87203
        ymin = 47.14912
        xmax = -52.60
        ymax = 47.33394
        if coordinates.geom_type.iloc[0] == 'Point':
            x = coordinates.x.iloc[0]
            y = coordinates.y.iloc[0]
        else:
            print(coordinates)
            return False
            raise ValueError("The coordinates are not of type 'Point' and cannot access .x/.y")
        x_within = (x > xmin) and (x < xmax)
        y_within = (y > ymin) and (y < ymax)
        if (not x_within) or (not y_within):
            print(y_within)
            print('puffin is not within the boundaries')
            exit(1)
            return False
        else:
            return True


    def get_line(self, point_in_line, multiline):
        coords = point_in_line.get_coordinates()
        x = coords.x.iloc[0]
        y = coords.y.iloc[0]
        minx = x - 0.00001
        miny = y - 0.00001
        maxx = x + 0.00001
        maxy = y + 0.00001
        line = multiline.clip((minx, miny, maxx, maxy))
        return line


    def compute_angle(self, line_1, line_2, direction):
        """
        Returns the smallest angle
        """
        l1 = line_1.get_coordinates()
        l2 = line_2.get_coordinates()
        v = np.array([l1.iloc[1].x-l1.iloc[0].x, l1.iloc[1].y-l1.iloc[0].y])
        u = np.array([l2.iloc[1].x-l2.iloc[0].x, l2.iloc[1].y-l2.iloc[0].y])
        unit_v = self.normalize_vector(v)
        unit_u = self.normalize_vector(u)
        angle = np.arccos(np.clip(np.dot(unit_v, unit_u), -1.0, 1.0))
        if direction > np.pi:
            return np.pi - angle
        else:
            return angle


    def normalize_vector(self, v):
        """ 
        Normalizes the given vector and returns it

        Parameters:
            - v (numpy.array): vector v to be normalized

        Returns:
            - u
        
        """
        u = v / np.linalg.norm(v)
        return u
