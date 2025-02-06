import geopandas as gpd
import matplotlib.pyplot as plt
from simulation import simulation


def main():
    crs = int(input('Enter crs to be used: '))
    puffins = int(input('Enter number of puffins: '))
    print()

    # Gull and Great islands, and Witless Bay coastline polyline geodataframe
    islands_polyline = gpd.read_file('/Users/flavio/Documents/puffin_simulation/files/2025-01-28_Puffin_colonies_shapefiles_FLT/islands_polyline.shp').to_crs(epsg=crs)
    witless_bay_polyline = gpd.read_file('/Users/flavio/Documents/puffin_simulation/files/coastline/coastline.shp').to_crs(epsg=crs)

    # Witless Bay polygon geodataframe
    witless_bay_polygon = gpd.read_file('/Users/flavio/Documents/puffin_simulation/files/Grid_puffmap_2km/Grid_puffmap_2km.shp').to_crs(epsg=crs)

    # Get Gull Island, Great Island and Witless Bay geoseries
    gull_polyline = islands_polyline[islands_polyline.Island=='gull']['geometry']
    great_polyline = islands_polyline[islands_polyline.Island=='great']['geometry']

    # run simulation
    sim_gull_gdf = simulation(gull_polyline, witless_bay_polyline, puffins, crs)
    sim_great_gdf = simulation(great_polyline, witless_bay_polyline, puffins, crs)

    # print results
    fig, ax = plt.subplots(figsize=(8, 8))
    sim_gull_gdf.plot(ax=ax, color='blue')
    sim_great_gdf.plot(ax=ax, color='green')
    witless_bay_polyline.plot(ax=ax, color='red', lw=0.2)
    plt.show()

if __name__ == '__main__':
    main()
