import geopandas as gpd
import numpy as np
import shapely
from shapely import Point, LineString
from puffin import Puffin


'''
This functions starts the simulation. The simulation starts by generating a puffin within a given
place, and then makes the puffin move until it reaches their destination.

Args:
    - start (geopandas.GeoSeries): GeoSeries object containing the geometry of the place where puffins spawn
    - destination (geopandas.GeoSeries): GeoSeries object containing the geometry of the destination
    - num_puffins (int): number of puffins to be used in the simulation
    - crs (int): coordinate reference system to be used in the simulation

Returns:
    - gdf (geopandas.GeoDataFrame): GeoDataFrame containing data about the puffins

'''
def simulation(start, destination, num_puffins, crs):

    # spawn puffins and store them in a array
    puffins = [spawn_puffin(start, crs) for i in range(num_puffins)]

    # initialized lists to be used to create a geodataframe
    position = []
    origin = []
    geometry = []       # shapely.Point list

    # start simulation
    in_progress = True

    while in_progress:

        # check if there are no puffins left in the list
        if len(puffins) == 0:
            in_progress = False
            continue

        puffin = puffins.pop()
        direction = puffin.get_direction()
        puffin_vector = puffin.get_line()

        # store puffin's initial position
        origin.append(puffin.get_position())

        # puffin just started moving
        is_moving = True

        while is_moving:
            dx = np.cos(direction) * 0.005              # change in x and y
            dy = np.sin(direction) * 0.005        
            puffin_vector = puffin_vector.translate(dx, dy)
            puffin.set_steps(puffin.get_steps()-1)

            # check if puffin has no steps left
            if puffin.get_steps() < 1:
                is_moving = False

            # check if puffin has reached their destination
            elif puffin_vector.intersects(destination).any():

                # update attributes of puffin
                intersection_point = puffin_vector.intersection(destination).iloc[0]
                new_coordinates = shapely.get_coordinates(intersection_point)[0]
                new_x = new_coordinates[0]
                new_y = new_coordinates[1]
                puffin.set_position((new_x, new_y))
                puffin.set_line(puffin_vector)

                # stop puffin
                is_moving = False

            if not is_moving:
                # add data to the respective lists
                position.append(puffin.get_position())
                geometry.append(Point(puffin.get_position()[0], puffin.get_position()[1]))

    # geodataframe                
    d = {'position': position, 'origin': origin, 'geometry': geometry}
    gdf = gpd.GeoDataFrame(d, crs=crs)
    
    return gdf

'''
This function creates an instance of the Puffin class.

Args:
    spawn (geopandas.GeoSeries): the geometry of the place where puffins spawn
    crs (int): coordinate reference system to be used

Returns:
    puffin (Puffin): instance of class Puffin
'''
def spawn_puffin(spawn, crs):

    steps = 100

    # get a random point within the geometry of the spawn
    origin = spawn.sample_points(size=1)
    
    # get a random direciton for the puffin
    direction = generate_angle((5*np.pi)/6, (7*np.pi)/6)

    # generate a vector for the puffin
    coords = origin.get_coordinates().iloc[0]
    x = coords['x']    
    y = coords['y']
    dx = x + np.cos(direction) * 0.005
    dy = y + np.sin(direction) * 0.005
    line = LineString([Point(dx, dy), Point(x, y)])

    # convert line to a geoseries object
    line = gpd.geopandas.GeoSeries(line, crs=crs)

    # initialize instance of puffin class
    puffin = Puffin((x, y), line, direction, steps, crs)

    return puffin

'''
This function takes an interval [min, max) and returns a angle in radians within such interval.

Args:
    - min, max (float): lowerbound and upperbound of the interval (both min and max must be in radians)

Returns:
    - angle (float): angle in [min, max)
'''
def generate_angle(min, max):

    angle = np.random.uniform(min, max)

    return angle
