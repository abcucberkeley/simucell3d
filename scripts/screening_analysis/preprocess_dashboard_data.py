import os
from sys import argv
import numpy as np
import pandas as pd
import json
from ast import literal_eval
from typing import Dict, Iterable, List, Literal, Tuple


"""
Utility functions to prepare data for plotting the dashboard.
"""

#-------------------------------------------------------------------------------------------------------
def clean_dashboard_files(
    path_to_dashboard_files_dir: str
) -> None:
    """
    Delete useless .lock files in the 'dashboard_files' folder and sort the files according to the
    simulation id.

    Parameters:
    -----------

    path_to_dashboard_files_dir: (str)
        The path to the directory where dashboard data are stored.        
    """

    # Remove .lock files
    for file in os.listdir(path_to_dashboard_files_dir):
        if file.endswith(".lock"):
            os.remove(os.path.join(path_to_dashboard_files_dir, file))
    
    # Sort files wrt sim_id
    fnames = os.listdir(path_to_dashboard_files_dir)
    for fname in fnames:
        data = [] # to temporary store the read content from the file
        
        if fname.endswith(".txt"):
            # read content
            with open(os.path.join(path_to_dashboard_files_dir, fname), "r") as f:
                for line in f:
                    data.append(line.strip().split(" "))
            # sort content
            sorted_data = sorted(data, key=lambda x: int(x[0]))
            # write sorted content
            with open(os.path.join(path_to_dashboard_files_dir, fname), "w") as fout:
                for idx, value in sorted_data:
                    fout.write(f"{idx} {value}\n")

        elif fname.endswith(".json"):
            # read content
            with open(os.path.join(path_to_dashboard_files_dir, fname), "r") as f:
                data = json.load(f)
            # sort content
            sorted_data = {key: data[key] for key in sorted(data.keys(), key=lambda x: int(x))}
            # write sorted content
            with open(os.path.join(path_to_dashboard_files_dir, fname), "w") as fout:
                json.dump(obj=sorted_data, fp=fout, indent=2)
#-------------------------------------------------------------------------------------------------------



#-------------------------------------------------------------------------------------------------------
def _map_parameters_to_table_names(
    param_names: Iterable[str]
) -> str:
    '''
    Given the name of a parameter, returns the actual name in the parameters table.
    Example: 'surface_tension' --> 'epi_api_surf_tens'

    Parameters:
    -----------
    param_names: (Iterable[str])
        The names of the parameters to show in the dashboard.
    
    Returns:
    --------
    table_names: (Iterable[str])
        The correspondent names in the parameters table.
    '''

    with open("./misc/param2table.json", "r") as file:
        param2table = json.load(file)

    return [param2table[param_name] for param_name in param_names]
#-------------------------------------------------------------------------------------------------------



#-------------------------------------------------------------------------------------------------------
def map_simulation_ids_to_parameters(
        screening_table_path: str,
        screened_parameters: Iterable[str],
        # save_dir: Optional[str] = None
) -> Dict[str, List[any]]:
    '''
    Extract from the screening parameter table the simulation ids and associate them to the 
    correspondent values of the screened parameters.

    Parameters:
    -----------
    screening_table_path: (str)
        Path to the .csv file used for setting the parameter screening.
    
    screened_parameters: Iterable[str]
        An array containing the names of the screened parameters.
        
    save_dir: (Optional[str], default=None)
        The directory in which the dictionary produced (.json file) by this function is saved.

    Returns:
    --------
    params_dict: (Dict[str, List[any]])
        A dictionary that associates 'sim_id' and names of screened parameters to their values.
    '''

    # Load csv file
    screen_table = pd.read_csv(screening_table_path)

    params_dict = {}

    #get simulation ids
    params_dict['sim_id'] = screen_table['sim_id'].to_numpy() 

    #get screened parameters columns
    table_names = _map_parameters_to_table_names(screened_parameters)
    for param, column in zip(screened_parameters, table_names):
        params_dict[param] = screen_table[column].to_numpy()

    return params_dict
#-------------------------------------------------------------------------------------------------------



#-------------------------------------------------------------------------------------------------------
def collect_indicators_data(
        indicators: List[Literal['delta_sphericity', 'interpenetration', 'volume_loss', 'num_iterations', 'IoU', 'IoU_derivative']],
        data_path: str,
) -> Dict[str, List[any]]:
    '''
    Collect previously computed values for indicators from the different simulation runs and store them in a 
    file as a dictionary indexed by simulation ids.

    Parameters:
    -----------
    indicators: Iterable[Literal['delta_sphericity', 'interpenetration', 'volume_loss', 'has_crashed']] DERIVATIVE OF SHAPE CHANGE
        A string defining an indicator name to include in the screening summary.
    
    data_path: (str)
        The path to the directory where the data files are stored.

    Returns:
    --------
    output_dict: (Dict[str, List[any]])
        A dictionary whose keys are 'sim_id' and the indicators names and the values are the correspondent data.
    '''

    # Clean and sort the dashboard files
    clean_dashboard_files(data_path)

    # Collect data from files
    out_dict = {}
    prev_sim_ids = []
    for indicator in indicators:
        assert indicator in ['delta_sphericity', 'interpenetration', 'volume_loss', 'num_iterations', 'IoU', 'IoU_derivative'], f"""\
            The selected indicator `{indicator}` is not defined.
            Please chose one among ['delta_sphericity', 'interpenetration', 'volume_loss', 'num_iterations', 'IoU', 'IoU_derivative'].
            """

        ext = 'json' if indicator in ['delta_sphericity', 'volume_loss', 'num_iterations'] else 'txt'
        data_file_path = os.path.join(data_path, f'{indicator}_output.{ext}')

        with open(data_file_path, 'r') as file:
            if ext == 'json':
                data_dict = json.load(file)
            elif ext == 'txt':
                data_dict = {}
                for line in file.readlines():
                    line = line.strip()
                    key, value = line.split(" ") 
                    data_dict[key] = literal_eval(value)

        # Check if sim_ids are consistent with previous file
        curr_sim_ids = sorted(list(data_dict.keys()), key=lambda x: int(x))
        if not prev_sim_ids:
            prev_sim_ids = curr_sim_ids
        else:
            assert curr_sim_ids == prev_sim_ids, "Simulation ids are not the same over different files" +\
                f"prev: {prev_sim_ids}, curr: {curr_sim_ids}"
            prev_sim_ids = curr_sim_ids

        data_dict = dict(sorted(data_dict.items(), key=lambda x: int(x[0])))
        data_lst = [val for val in data_dict.values()]
        out_dict[indicator] = data_lst

    out_dict['sim_id'] = [int(sim_id) for sim_id in prev_sim_ids]

    return out_dict
#-------------------------------------------------------------------------------------------------------



#-------------------------------------------------------------------------------------------------------
def merge_dataframes(
        path_to_screening_table: str,
        path_to_data_dir: str,
        param_names: List[str],
        ind_names: List[Literal['delta_sphericity, interpenetration', 'volume_loss', 'num_iterations', 'IoU', 'IoU_derivative']]
) -> pd.DataFrame:
    '''
    Load parameters from screening table and previously computed indicator data from files.
    Merge them in a dataframe and return it.

    Parameters:
    -----------
    path_to_screening_table: (str)
        Path to the .csv file used for setting the parameter screening.
    
    path_to_data_dir: (str)
        The path to the directory where the data files are stored.
    
    param_names: (List[str])
        An array containing the names of the screened parameters.
    
    ind_names: (List[Literal['delta_sphericity, interpenetration', 'volume_loss', 'num_iterations', 'IoU', 'IoU_derivative']])
        A string defining an indicator name to include in the screening summary.

    Returns:
    --------
    res_df: (pd.DataFrame)
        A dataframe that gathers parameters from screening table and previously computed indicator data from files.
    '''
    
    params_dict = map_simulation_ids_to_parameters(
        screening_table_path=path_to_screening_table,
        screened_parameters=param_names
    )

    indicators_dict = collect_indicators_data(
        indicators=ind_names,
        data_path=path_to_data_dir
    )

    # Get the common simulation ids
    common_sim_ids_mask = np.isin(params_dict['sim_id'], indicators_dict['sim_id'])
    for key in params_dict.keys():
        params_dict[key] = [val for val, mask in zip(params_dict[key], common_sim_ids_mask) if mask]

    assert np.all(params_dict['sim_id'] == indicators_dict['sim_id']), f"""\
        No correspondence in simulation ids! {params_dict['sim_id']} {indicators_dict['sim_id']}"""

    params_df = pd.DataFrame(params_dict)
    indicators_df = pd.DataFrame(indicators_dict)
    res_df = pd.merge(params_df, indicators_df, on='sim_id', how='outer')

    return res_df
#-------------------------------------------------------------------------------------------------------



#-------------------------------------------------------------------------------------------------------
def _split_indicator_cols(
    df: pd.DataFrame
) -> pd.DataFrame:
    '''
    Given the dataframe of the data to represent in the dashboard, parse its column
    splitting the columns whose elements are list/tuples (e.g., mean + std or num + bool).

    Parameters:
    -----------
    df: (pd.DataFrame)
        The dataframe storing the data to report in the dashboard.

    Returns:
    --------
    out_df: (pd.DataFrame)
        A copy of the input dataframe with splitted columns.
    ''' 

    out_df = df.copy()

    for column in out_df.columns:
        if column in ['delta_sphericity', 'volume_loss']:
            means = out_df[column].apply(lambda val: val[0] if isinstance(val, list) and len(val) > 0 else np.nan)
            stds = out_df[column].apply(lambda val: val[1] if isinstance(val, list) and len(val) > 0 else np.nan)
            out_df.drop(columns=[column], inplace=True)
            out_df[column + '_mean'], out_df[column + '_std'] = means, stds
        elif column == 'num_iterations':
            vals = out_df[column].apply(lambda val: val[0] if isinstance(val, list) and len(val) > 0 else np.nan)
            logics = out_df[column].apply(lambda val: val[1] if isinstance(val, list) and len(val) > 0 else np.nan)
            out_df.drop(columns=[column], inplace=True)
            out_df[column], out_df['has_crashed'] = vals, logics
    
    return out_df
#-------------------------------------------------------------------------------------------------------



#-------------------------------------------------------------------------------------------------------
def _find_failed_simulation_runs(
    df: pd.DataFrame,   
) -> pd.DataFrame:
    '''
    Given the dataframe of the data to represent in the dashboard, find the records corresponding
    to FAILED simulation runs. We have two types of those:
    - Runs that have crashed -> `num_iterations`
    - Runs that lose one cell after some iterations -> `None` in `volume_loss`/`sphericity`

    Parameters:
    -----------
    df: (pd.DataFrame)
        The dataframe storing the data to report in the dashboard.

    Returns:
    --------
    df: (pd.DataFrame)
        The input dataframe with an additional boolean column `has_failed`.
    '''
    
    failed_lst = []
    for _, row in df.iterrows():
        new_val = np.isnan(row['delta_sphericity_mean']) or row['has_crashed']
        failed_lst.append(new_val)
    
    df['has_failed'] = failed_lst
    return df
#-------------------------------------------------------------------------------------------------------



#-------------------------------------------------------------------------------------------------------
def _find_outliers(
    df: pd.DataFrame,   
) -> pd.DataFrame:
    '''
    Given the dataframe of the data to represent in the dashboard, find the outlying records 
    for volume loss and delta sphericity.

    Parameters:
    -----------
    df: (pd.DataFrame)
        The dataframe storing the data to report in the dashboard.

    Returns:
    --------
    df: (pd.DataFrame)
        The input dataframe with an additional boolean column `is_outlier`.
    '''
    # Find `sphericity` and `volume_loss` outliers
    outliers = np.zeros(len(df), dtype=bool)
    for column in ['delta_sphericity_mean', 'volume_loss_mean']:
        data = abs(df[column].values)
        upper_bound = np.quantile(data, 0.99)
        outliers = np.logical_or(outliers, (data > upper_bound))
                
    df['is_outlier'] = outliers.astype(bool)
    return df
#-------------------------------------------------------------------------------------------------------



#-------------------------------------------------------------------------------------------------------
def preprocess_data(
    df: pd.DataFrame,
    parameters: Iterable[str],
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    '''
    Prepare dataframe for the dashboard. The cleaning involves the followig steps:
    - Split the indicator columns with multiple values,
    - Add column to mark failed simulation runs
    - Return a clean dataframe with only successful runs (no NaNs) and one containing grid
      of parameters corresponding to failed runs.
    
    Parameters:
    -----------
    df: (pd.DataFrame)
        The dataframe storing the data to report in the dashboard.

    parameters: (Iterable[str])
        The list of screening parameters.

    Returns:
    --------
    clean_df: (pd.DataFrame)
        A copy of the input dataframe without data associated to failed simulation runs.

    failed_df: (pd.DataFrame)
        A dataframe with only parameter columns associated to failed simulation runs.
    
    outliers_df: (pd.DataFrame)
        A dataframe with only parameter columns associated to outlying values forr indicators.

    interpenetration_df: (pd.DataFrame)
        A dataframe with only parameter columns associated to runs in which interpenetration was spotted.
    '''

    out_df = df.copy()
    out_df = _split_indicator_cols(out_df)

    # Mark failed and outlying simulation runs
    out_df = _find_failed_simulation_runs(out_df)
    out_df = _find_outliers(out_df)

    # Isolate dataframes of parameters associated to failed, outlying, and interpenetrating runs
    failed_df = out_df.loc[out_df['has_failed'], parameters].copy()
    outliers_df = out_df.loc[out_df['is_outlier'], parameters].copy()
    interpenetration_df = out_df.loc[out_df['interpenetration'], parameters].copy()
    
    # Remove from original data all the 'dirty' runs
    to_keep = np.logical_and(
        np.logical_and(~out_df['is_outlier'].values, ~out_df['interpenetration'].values),
        ~out_df['has_failed'].values
    )
    clean_df = out_df[to_keep]


    return clean_df, failed_df, outliers_df, interpenetration_df
#-------------------------------------------------------------------------------------------------------