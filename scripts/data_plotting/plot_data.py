import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sys import argv
from os import path, _exit, listdir
from tqdm import tqdm
import warnings

warnings.filterwarnings("ignore", category=RuntimeWarning)

"""
    This script plots the data from the simulation_statistics.csv file. The plots will be saved 
    in the same folder as the simulation_statistics.csv file.
"""



def get_avg_and_std_of_data(grouped_data, column_name):
    """
        This function calculates the average and standard deviation of a given data array.

        Parameters:
        -----------
        
        grouped_data (pd.DataFrameGroupBy()): 
            The data grouped by iteration.
        
        column_name (str):
            The name of the column to get the data from.

        Returns:
        --------

        avg (np.array()):
            The average of the column for each iteration

        std (np.array()):
            The standard deviation of the column for each iteration
    """


    avg = grouped_data[column_name].mean().to_numpy()
    std = grouped_data[column_name].std().to_numpy()
    return avg, std





def get_deviation_from_target_value(grouped_data, column_value, column_target_value):
    """
        This function calculates the deviation of the cell from their target areas or target volumes.

        Parameters:
        -----------
        
        grouped_data (pd.DataFrameGroupBy()): 
            The data grouped by iteration.
        
        column_name (str):
            The name of the column to get the data from (volume or area)

        column_target_value (str):
            The name of the column to get the target value from (target_volume or target_area)

        Returns:
        --------

        avg (np.array()):
            The average of the column for each iteration

        std (np.array()):
            The standard deviation of the column for each iteration
    """

    #Store the average and std of the deviation from the target value in these lists
    avg_deviation_from_target_value = []
    std_deviation_from_target_value = []

    for name, group in grouped_data:

        #Get the deviation from the target value
        deviation_from_target_value = np.absolute(group[column_value].to_numpy() - group[column_target_value].to_numpy()) / group[column_target_value].to_numpy()

        #Calculate the average and std of the deviation from the target value
        avg_deviation_from_target_value.append(np.mean(deviation_from_target_value))
        std_deviation_from_target_value.append(np.std(deviation_from_target_value))


    return np.array(avg_deviation_from_target_value), np.array(std_deviation_from_target_value)


def plot_stats(
        path_to_stats_df: str,
        save_dir: str
):  
    '''
    Extract simulation statistics from simulation runs dataframes and plot their evolution.

    Parameters:
    -----------
        path_to_stats_df: (str)
            The path to the 'simulation_statistics_{idx}.csv file.

        save_dir: (str)
            The path to the directory where to save the plots.
    '''

    # Read in the data
    sim_data = pd.read_csv(path.join(path_to_stats_df))

    #Get the start and end time of the simulation
    start_time = sim_data['simulation_time'].min()
    end_time =   sim_data['simulation_time'].max()

    #Group the data by iteration
    sim_data_grouped = sim_data.groupby('iteration')


    #Get the number of cells as a function of time
    nb_cells_ar = sim_data_grouped.size().to_numpy()

    #Get the different arrays of data
    avg_vol, std_vol                           = get_avg_and_std_of_data(sim_data_grouped, 'volume')
    avg_dev_target_vol, std_dev_target_vol     = get_deviation_from_target_value(sim_data_grouped, 'volume', 'target_volume')
    avg_area, std_area                         = get_avg_and_std_of_data(sim_data_grouped, 'area')
    avg_dev_target_area, std_dev_target_area   = get_deviation_from_target_value(sim_data_grouped, 'area', 'target_area')
    avg_pressure, std_pressure                 = get_avg_and_std_of_data(sim_data_grouped, 'pressure')
    

    avg_kinetic_energy, std_kinetic_energy                          = get_avg_and_std_of_data(sim_data_grouped, 'kinetic_energy')
    avg_surface_tension_energy, std_surface_tension_energy          = get_avg_and_std_of_data(sim_data_grouped, 'surface_tension_energy')
    avg_membrane_elasticity_energy, std_membrane_elasticity_energy  = get_avg_and_std_of_data(sim_data_grouped, 'membrane_elasticity_energy')
    avg_bending_energy, std_bending_energy                          = get_avg_and_std_of_data(sim_data_grouped, 'bending_energy')
    avg_pressure_energy, std_pressure_energy                        = get_avg_and_std_of_data(sim_data_grouped, 'pressure_energy')
    avg_total_energy, std_total_energy                              = get_avg_and_std_of_data(sim_data_grouped, 'total_potential_energy')

    #Get the simulation time at each iteration
    simulation_time_ar = sim_data_grouped['simulation_time'].min().to_numpy()

    #Plot the data related to the grometries and to the energies in 2 different figures
    fig_1 = plt.figure(figsize=(8, 6))

    ax_1 = fig_1.add_subplot(2, 3, 1)
    ax_1.plot(simulation_time_ar, nb_cells_ar)
    ax_1.set_xlabel('Simulation time (s)')
    ax_1.set_ylabel('Number of cells')
    ax_1.set_xlim(start_time, end_time)

    #Plot the volume of the cell
    ax_2 = fig_1.add_subplot(2, 3, 2)
    ax_2.fill_between(simulation_time_ar, avg_vol-std_vol, avg_vol+std_vol, alpha = 0.5)
    ax_2.plot(simulation_time_ar, avg_vol)
    ax_2.set_xlabel('Simulation time (s)')
    ax_2.set_ylabel('Cell volume (m^3)')
    ax_2.set_xlim(start_time, end_time)

    #Plot the deviation from the target volume of the cell
    ax_3 = fig_1.add_subplot(2, 3, 3)
    ax_3.fill_between(simulation_time_ar, avg_dev_target_vol-std_dev_target_vol, avg_dev_target_vol+std_dev_target_vol, alpha = 0.5)
    ax_3.plot(simulation_time_ar, avg_dev_target_vol)
    ax_3.set_xlabel('Simulation time (s)')
    ax_3.set_ylabel('Deviation from target\nvolume |(v - v_0)| / v_0')
    ax_3.set_xlim(start_time, end_time)

    #Plot the area of the cell
    ax_4 = fig_1.add_subplot(2, 3, 4)
    ax_4.fill_between(simulation_time_ar, avg_area-std_area, avg_area+std_area, alpha = 0.5)
    ax_4.plot(simulation_time_ar, avg_area)
    ax_4.set_xlabel('Simulation time (s)')
    ax_4.set_ylabel('Cell area (m^2)')
    ax_4.set_xlim(start_time, end_time)

    #Plot the deviation from the target area of the cell
    ax_5 = fig_1.add_subplot(2, 3, 5)
    ax_5.fill_between(simulation_time_ar, avg_dev_target_area-std_dev_target_area, avg_dev_target_area+std_dev_target_area, alpha = 0.5)
    ax_5.plot(simulation_time_ar, avg_dev_target_area)
    ax_5.set_xlabel('Simulation time (s)')
    ax_5.set_ylabel('Deviation from target\narea |(a - a_0)| / a_0')
    ax_5.set_xlim(start_time, end_time)

    #Plot the pressure of the cell
    ax_6 = fig_1.add_subplot(2, 3, 6)
    ax_6.fill_between(simulation_time_ar, avg_pressure-std_pressure, avg_pressure+std_pressure, alpha = 0.5)
    ax_6.plot(simulation_time_ar, avg_pressure)
    ax_6.set_xlabel('Simulation time (s)')
    ax_6.set_ylabel('Cell pressure (Pa)')
    ax_6.set_xlim(start_time, end_time)

    fig_1.tight_layout()
    fig_1.savefig(path.join(save_dir, "sim_geometry_data.png"))
    plt.clf()

    #Plot the different energies of the cells
    fig_2 = plt.figure(figsize=(8, 6))
    
    ax_7 = fig_2.add_subplot(2, 3, 1)
    ax_7.fill_between(simulation_time_ar, avg_kinetic_energy-std_kinetic_energy, avg_kinetic_energy+std_kinetic_energy, alpha = 0.5)
    ax_7.plot(simulation_time_ar, avg_kinetic_energy)
    ax_7.set_xlabel('Simulation time (s)')
    ax_7.set_ylabel('Kinetic energy (J)')
    ax_7.set_yscale('log')
    ax_7.set_xlim(start_time, end_time)

    ax_8 = fig_2.add_subplot(2, 3, 2)
    ax_8.fill_between(simulation_time_ar, avg_surface_tension_energy-std_surface_tension_energy, avg_surface_tension_energy+std_surface_tension_energy, alpha = 0.5)
    ax_8.plot(simulation_time_ar, avg_surface_tension_energy)
    ax_8.set_xlabel('Simulation time (s)')
    ax_8.set_ylabel('Surface tension energy (J)')
    ax_8.set_yscale('log')
    ax_8.set_xlim(start_time, end_time)


    ax_9 = fig_2.add_subplot(2, 3, 3)
    ax_9.fill_between(simulation_time_ar, avg_membrane_elasticity_energy-std_membrane_elasticity_energy, avg_membrane_elasticity_energy+std_membrane_elasticity_energy, alpha = 0.5)
    ax_9.plot(simulation_time_ar, avg_membrane_elasticity_energy)
    ax_9.set_xlabel('Simulation time (s)')
    ax_9.set_ylabel('Membrane elasticity energy (J)')
    # ax_9.set_yscale('log')
    ax_9.set_xlim(start_time, end_time)

    ax_10 = fig_2.add_subplot(2, 3, 4)
    ax_10.fill_between(simulation_time_ar, avg_pressure_energy-std_pressure_energy, avg_pressure_energy+std_pressure_energy, alpha = 0.5)
    ax_10.plot(simulation_time_ar, avg_pressure_energy)
    ax_10.set_xlabel('Simulation time (s)')
    ax_10.set_ylabel('Pressure energy (J)')
    ax_10.set_yscale('log')
    ax_10.set_xlim(start_time, end_time)

    ax_11 = fig_2.add_subplot(2, 3, 5)
    ax_11.fill_between(simulation_time_ar, avg_total_energy-std_total_energy, avg_total_energy+std_total_energy, alpha = 0.5)
    ax_11.plot(simulation_time_ar, avg_total_energy)
    ax_11.set_xlabel('Simulation time (s)')
    ax_11.set_ylabel('Total energy (J)')
    ax_11.set_yscale('log')
    ax_11.set_xlim(start_time, end_time)

    fig_2.tight_layout()
    fig_2.savefig(path.join(save_dir, "sim_energy_data.png"))
    plt.clf()


if __name__ == "__main__":

    import sys

    assert sys.argv.__len__() == 2, r"Script usage: python3 plot_data path/to/simulation_output_folder"


    
    sim_stats_folder = sys.argv[1]
    plot_stats(
        path_to_stats_df=path.join(sim_stats_folder, 'simulation_statistics.csv'),
        save_dir = sim_stats_folder
    )







    
    