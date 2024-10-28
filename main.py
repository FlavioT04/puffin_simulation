import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import simulation
import numpy as np


def main():
    crs = int(input('Enter the crs to be used: '))
    number_of_puffins = int(input('Enter number of puffins: '))
    print()

    # Base map
    provinces = gpd.read_file('files/provinces/lpr_000b21a_e.shp').to_crs(epsg=crs)
    newfoundland = provinces[provinces.PRNAME == 'Newfoundland and Labrador / Terre-Neuve-et-Labrador']
    newfoundland = gpd.clip(newfoundland, (-52.8395, 47.24035, -52.7624, 47.29752))

    # Witless Bay map
    witless_bay = gpd.GeoDataFrame(geometry=[newfoundland.geometry.iloc[0].geoms[5]], crs=crs)

    # Gull Island map
    map = gpd.read_file('files/census_divisions/lcsd000b21a_e.shp')
    map = map[map.CSDNAME == 'Witless Bay'].to_crs(epsg=crs)
    gull_island = gpd.GeoDataFrame(geometry=[map.geometry.iloc[0].geoms[0]], crs=crs)
 
    # Run simulation
    sim_1 = simulation.Simulation(gull_island, witless_bay, number_of_puffins, crs)
    result = sim_1.start()

    # Puffins GeoDataFrame
    puffin_data = {'Puffin': result}
    puffin_gdf = gpd.GeoDataFrame(puffin_data, geometry=result, crs=crs)

    # Spatial join
    map_with_puffins = gpd.sjoin(newfoundland, puffin_gdf, how='inner', predicate='contains')
    base = map_with_puffins.plot(figsize=(12, 8), color='#D3D3D3', edgecolor='black', linewidth=0.5)

    # Plotting
    ax = puffin_gdf.plot(ax=base, marker='o', color='black', markersize=5)
    ax.set_xlim([-52.83935, -52.74986])
    ax.set_ylim([47.24039, 47.29664])
    plt.show()


if __name__ == '__main__':
    main()