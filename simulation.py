import geopandas as gpd
import numpy as np
from shapely import Point, LineString
from puffin import Puffin
import matplotlib.pyplot as plt


def simulation(start, destination, num_puffins, crs):

    # spawn puffins and store them in a array
    puffins = [spawn_puffin(start, crs) for i in range(num_puffins)]

    # initialized lists to be used to create a geodataframe
    position = []
    origin = []
    geometry = []

    # start simulation
    in_progress = True

    while in_progress:
        total_dx = 0        # the total change in x and y
        total_dy = 0

        # check if there are no puffins left in the list
        if len(puffins) == 0:
            in_progress = False
            continue

        puffin = puffins.pop()
        direction = puffin.get_direction()
        path = puffin.get_line()

        # store puffin's initial position
        origin.append(puffin.get_position())

        # puffin just started moving
        is_moving = True

        while is_moving:
            dx = np.cos(direction) * 0.0006              # change in x and y
            dy = np.sin(direction) * 0.0006
            total_dx += dx
            total_dy += dy           
            path = path.translate(dx, dy)

            # check if puffin has reached its destination
            if path.intersects(destination).any():
                is_moving = False
                print('puffin reached destination\n')

                # update attributes of puffin
                x = puffin.get_position()[0]
                y = puffin.get_position()[1]
                puffin.set_position((x + total_dx, y + total_dy))
                puffin.set_line(path)

                # add data to the respective lists
                position.append(puffin.get_position())
                geometry.append(puffin.get_line())

                continue
    
    return

'''
This function gets everything necessary to create an instance of the puffin class.

Args:
    geometry (geopandas.GeoSeries): geometry used to retrive a point within itself
    crs (int): coordinate reference system to be used

Returns:
    puffin (Puffin): instance of class Puffin
'''
def spawn_puffin(geometry, crs):

    # get random point within geometry
    origin = geometry.sample_points(size=1)
    
    # get random direciton for puffin
    direction = compute_random_direction((5*np.pi)/6, (7*np.pi)/6)

    # create line to represent path of the puffin
    coords = origin.get_coordinates()
    x = coords['x']    
    y = coords['y']
    dx = x + np.cos(direction) * 0.005
    dy = y + np.sin(direction) * 0.005
    line = LineString([Point(dx, dy), Point(x, y)])

    # convert line to a geoseries object
    line = gpd.geopandas.GeoSeries(line, crs=crs)

    # initialize instance of puffin class
    puffin = Puffin((x, y), line, direction, crs)

    return puffin

'''
This function computes a floating point number representing a direction, in randins, within the interval [min, max) taken
by a puffin. The direction is determined randomly.

Args:
    - min, max (float): lowerbound and upperbound of the interval [min, max)

Returns:
    - direction (float): direction taken by puffin
'''
def compute_random_direction(min, max):

    direction = np.random.uniform(min, max)

    return direction
