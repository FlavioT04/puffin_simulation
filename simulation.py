import geopandas as gpd
import numpy as np
import shapely
import matplotlib.pyplot as plt
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

    # files needed
    witless_bay_grid = gpd.read_file('files/Grid_puffmap_2km/Grid_puffmap_2km.shp').to_crs(epsg=crs)

    # settings
    b_height = 6
    b_angle = np.pi/4

    # modes
    bouncing_angle = False
    bouncing_wall = True

    # spawn puffins and store them in a array
    puffins = [spawn_puffin(start, crs) for i in range(num_puffins)]

    # initialized lists to be used to create a geodataframe
    position = []
    origin = []
    geometry = []       # shapely.Point list
    vector_geometry = []

    # start simulation
    in_progress = True

    while in_progress:

        # check if there are no puffins left in the list
        if len(puffins) == 0:
            in_progress = False
            continue

        puffin = puffins.pop()
        puffin_vector = puffin.get_vector()

        # store puffin's initial position
        origin.append(puffin.get_position())

        # puffin just started moving
        is_moving = True

        while is_moving:
            
            # move puffin
            puffin_vector = puffin.move()

            # check if puffin has no steps left
            if puffin.get_steps() < 1:
                print('-------------------------------')
                print('puffin killed')
                print('-------------------------------')
                is_moving = False

            # check if puffin has reached their destination
            elif puffin_vector.intersects(destination).any():

                # update vector
                puffin.set_vector(puffin_vector)

                # when bouncing_wall is True, puffin bounces depending on the elevation of the intersecting geometry
                if bouncing_wall == True:
                    
                    # get the elevation of the cell
                    cell_elevation = intersecting_cell(witless_bay_grid, puffin_vector)['MeanElevat']

                    # check conditions
                    if cell_elevation <= b_height:
                        print('------Puffin moves inland------')
                        print('cell elevation:', cell_elevation)
                        puffin.move()
                        while is_moving:
                            # get random integer from 0 to 100 (exclusive)
                            random_num = np.random.randint(0, 100)
                            print(random_num)        
                            if cell_elevation >= 0 and cell_elevation <= 2:
                                if random_num < 50:
                                    is_moving = False
                                else: 
                                    for i in range(20):
                                        puffin_vector = puffin.move()
                            if cell_elevation > 2 and cell_elevation <= 4:
                                if random_num < 75:
                                    is_moving = False
                                else: 
                                    for i in range(20):
                                        puffin_vector = puffin.move()
                            if cell_elevation > 4:
                                is_moving = False
                            if cell_elevation < 0:
                                print('elevation is negative')
                                is_moving = False
                            
                            if is_moving:
                                cell = intersecting_cell(witless_bay_grid, puffin_vector)
                                if cell['cell_ID'] == 'None':
                                    print('new cell elevation:', None)
                                    is_moving = False
                                else:
                                    cell_elevation = cell['MeanElevat']
                                print('new cell elevation:', cell_elevation)
                    else:
                        print('--------puffin bounces--------')
                        print('cell elevation:', cell_elevation)
        
                # stop puffin
                is_moving = False

            if not is_moving:
                if puffin.get_steps() > 0:
                    print('puffin stopped')
                # add data to the respective lists
                position.append(puffin.get_position())
                vector_geometry.append(puffin.get_vector().iloc[0])
                geometry.append(Point(puffin.get_position()[0], puffin.get_position()[1]))
                print()
        

    # geodataframe                
    d = {'position': position, 'origin': origin, 'geometry': geometry}
    gdf = gpd.GeoDataFrame(d, crs=crs)

    # vector geodataframe
    v_d = {'geometry': vector_geometry}
    v_gdf = gpd.GeoDataFrame(v_d, crs=crs)
    
    return gdf, v_gdf

'''
This function creates an instance of the Puffin class.

Args:
    spawn (geopandas.GeoSeries): MultiLineString geometry of the place where puffins spawn
    crs (int): coordinate reference system to be used

Returns:
    puffin (Puffin): instance of class Puffin
'''
def spawn_puffin(spawn, crs):

    # num of steps a puffin can take
    steps = 5000

    # get a random point within the geometry of the spawn
    origin = spawn.sample_points(size=1)

    # clip the line containing the origin
    intersecting_line = clip_geometry(spawn, origin).set_crs(crs)

    # generate a vector for the puffin
    coords = origin.get_coordinates().iloc[0]
    x = coords['x']    
    y = coords['y']
    dx = x + 0.00005 * np.cos(np.pi)
    dy = y
    vector = gpd.geopandas.GeoSeries(LineString([Point(dx, dy), Point(x, y)])).set_crs(epsg=crs)

    # get angle between the puffin's vector and the spawn
    angle = get_intersecting_angle(vector, intersecting_line)
    
    # rotate vector so that it is perpendicular to the LineString containing its spawn point
    alpha = np.abs(np.pi/2 - angle)
    if angle >= np.pi/2:
        vector = vector.rotate(angle=alpha, origin=origin.iloc[0], use_radians=True)
    else:
        vector = vector.rotate(angle=-alpha, origin=origin.iloc[0], use_radians=True)
        
    # rotate vector by random angle
    vector = rotate_by_random_angle(5*np.pi/6, 7*np.pi/6, vector, origin.iloc[0])

    # get the direction of the vector
    vector_coords = vector.get_coordinates()
    direction = np.arctan((vector_coords.iloc[1].y-vector_coords.iloc[0].y) / (vector_coords.iloc[1].x-vector_coords.iloc[0].x))
    if direction < np.pi:
        direction = np.pi + direction

    # Puffin
    puffin = Puffin((x, y), vector, direction, steps, crs)

    return puffin

'''
This function takes an interval [min, max) and returns a angle in radians within such interval.

Args:
    - min, max (float): lowerbound and upperbound of the interval (both min and max must be in radians)
    - vector (geopandas.GeoSeries): vector of the puffin
    - point (shapely.Point): origin of rotation

Returns:
    - rotated_vector (geopandas.Geoseries): rotated vector
'''
def rotate_by_random_angle(min, max, vector, point):
    
    if max < min:
        print('the maximum angle is smaller than the minimum angle')
        exit(1)

    # get distance between angles
    da = max - min

    # rotate by random angle
    random_angle = np.random.default_rng().uniform(-da, da)
    rotated_vector = vector.rotate(angle=random_angle, origin=point, use_radians=True)

    return rotated_vector

'''
This function computes the angle between two intersecting vectors

Args:
    - v, u (geopandas.GeoSeries): geoseries objects with LineStrings as their geometry

Returns:
    - angle (float): angle of intersection
'''
def get_intersecting_angle(v, u):
    
    l1 = v.get_coordinates()
    l2 = u.get_coordinates()

    # create numpy array
    vn = np.array([l1.iloc[1].x-l1.iloc[0].x, l1.iloc[1].y-l1.iloc[0].y])
    un = np.array([l2.iloc[1].x-l2.iloc[0].x, l2.iloc[1].y-l2.iloc[0].y])

    # compute magnitudes
    magv = np.sqrt(vn[0]**2 + vn[1]**2)
    magu = np.sqrt(un[0]**2 + un[1]**2)

    # compute angle
    angle = np.arccos(np.dot(vn, un) / (magv * magu))

    return angle

'''
This function clips a geometry

Args:
    - geometry (geopandas.GeoSeries): geometry to be clipped
    - point (geopandas.GeoSeries): point in geometry

Returns:
    - clipped (geopandas.GeoSeries): clipped geometry
'''
def clip_geometry(geometry, point):

    x = point.x.iloc[0]
    y = point.y.iloc[0]
    minx = x - 0.00001
    miny = y - 0.00001
    maxx = x + 0.00001
    maxy = y + 0.00001
    clipped = geometry.clip((minx, miny, maxx, maxy))

    return clipped

'''
This function looks for the cell containing the given puffin

Args:
    - grid (geopandas.GeoDataFrame): grid to be used in the search
    - puffin_vector (geopandas.GeoSeries): geometry of the puffin

Returns:
    - cell (geopandas.GeoSeries): cell containing the given puffin
'''
def intersecting_cell(grid, puffin_vector):
    searching = True
    cell = gpd.GeoDataFrame({'cell_ID': ['None'], 'geometry': Point(0, 0)}, crs=grid.crs).iloc[0]
    i = 0
    while searching and i < len(grid):
        cell_geometry = grid.iloc[i]['geometry']
        if cell_geometry.intersects(puffin_vector.iloc[0]):
            cell = grid.iloc[i]
            searching = False
        else:
            i += 1
    return cell        
