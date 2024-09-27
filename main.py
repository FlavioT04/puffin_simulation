from puffin_simulation import simulate_puffin_movement
from land_generator import LandGenerator

def main():
    gull_island = LandGenerator("Gull Island", 0, 0)
    witless_bay = LandGenerator("Witless Bay", -1, 0)
    print("source:", gull_island)
    simulate_puffin_movement(gull_island, witless_bay)

if __name__ == "__main__":
    main()