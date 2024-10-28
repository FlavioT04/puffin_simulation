import pandas as pd
import geopandas as gpd
import numpy as np


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
            puffin = self.spawn_puffin(self.source)    # GeoSeries Point object
            direction = self.random_angle((5*np.pi)/6, (7*np.pi/6))    # Angle within [150, 210)
            destination_reached = False
            while not destination_reached:
                puffin = self.move_puffin(puffin, direction)
                puffin_coordinates = puffin.get_coordinates()
                puffin_x = puffin_coordinates.x
                if (puffin_x < -53).any() == True:
                    break
                check = puffin.within(self.destination.geometry)
                if check.all() == False:
                    destination_reached = False
                else:
                    destination_reached = True
                    puffins_reached_destination +=1
                    self.puffins.append(puffin[0])
                    print(f'Puffin {puffins_reached_destination} reached destination')
                    print()
        print('Number of puffins:', len(self.puffins))
        print()
        return self.puffins


    def get_random_coordinate(self, space):
        """
        Returns a GeoSeries Point with coordinates within the geometry of the given space

        Parameters:
            - space (geopandas.GeoSeries): GeoSeries of the space

        Returns:
            - coordinate (shapely.Point): Point with coordinates x, y
        """
        random_point = space.sample_points(size=1)
        coordinate = random_point
        return coordinate


    def spawn_puffin(self, origin):
        """
        Spawns a puffin within the specified origin

        Parameters:
            - origin (geopandas.GeoDataFrame): GeoDataFrame containing geometry where puffin is to be spawned within

        Returns:
            - new_puffin (geopandas.GeoSeries): Point with coordinates (x, y) within the origin
        """
        new_puffin = gpd.GeoSeries(self.get_random_coordinate(origin), crs=self.crs)
        return new_puffin
    

    def move_puffin(self, puffin, direction):
        """
        Translates the geometry of a puffin

        Parameters:
            - puffin (GeoSeries Point): Point representing puffin
            - direction (float): angle of direction

        Returns:
            - translated_puffin (geopandas.GeoSeries): Point with translated geometry
        """
        dx = np.cos(direction) * 0.00006
        dy = np.sin(direction) * 0.00006
        translated_puffin = gpd.GeoSeries.translate(puffin, dx, dy)
        return translated_puffin

    
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
