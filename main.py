import sys
import pandas as pd
from puffin_simulation import simulate_puffin_movement
from land_generator import LandGenerator

def main(csv_filename, emitter_name, receiver_name):
    """
    Main function to initialize the puffin simulation. It reads the CSV file, extracts coordinates for the emitter and
    receiver, and initializes the simulation

    args:
        csv_filename (str): Path to the csv file containing the coordinates of the detectors
        emitter_name (str): Name of the land where the emitters are located
        receiver_name (str): Name of the land where the receivers are located

    The CSV file should contain the following columns:
        - x (dd): x-coordinates
        - y (dd): y-coordinates
        - role: Indicates whether the row represents an 'emitter' or 'receiver'
    """

    # Read CSV file into a pandas dataframe
    df = pd.read_csv(csv_filename)

    # List containing tuples (x, y) to represent the coordinates of the emitter
    emitter_coords = df[df['role'] == 'emitter']
    emitter_coords = [(row['x (dd)'], row['y (dd)']) for _, row in emitter_coords.iterrows()]

    # List containing tuples (x, y) to represent the coordinates of the receiver
    receiver_coords = df[df['role'] == 'receiver']
    receiver_coords = [(row['x (dd)'], row['y (dd)']) for _, row in receiver_coords.iterrows()]

    # Initialize simulation
    emitter_land = LandGenerator(emitter_name, 'emitter', emitter_coords)
    receiver_land = LandGenerator(receiver_name, 'receiver', receiver_coords)
    print(f"Starting simulation...\nSource: {emitter_land}\nDestination: {receiver_land}\n")
    simulate_puffin_movement(emitter_land, receiver_land) 

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Usage: python main.py <csv_filename>')
    else:
        csv_filename= sys.argv[1]
        emitter_name = sys.argv[2]
        receiver_name = sys.argv[3]
        main(csv_filename, emitter_name, receiver_name)
