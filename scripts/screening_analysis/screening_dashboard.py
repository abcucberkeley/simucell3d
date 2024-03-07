import os
import json
import numpy as np
import pandas as pd
import matplotlib
from matplotlib.colors import to_rgba
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, LogFormatter
import seaborn as sns
from itertools import combinations
from typing import Optional, Dict, Iterable, Literal, Tuple
from preprocess_dashboard_data import merge_dataframes, preprocess_data


"""
Given a parameter screening made by a bunch of simulation runs, create a dashboard to summarize the outcomes.
The dashboard is composed by plots reporting key summary statistics (e.g., volume contraction, 
sphericty at end of simulation) and indicators (e.g., simulation has crashed, interpenetration, ...).
"""


#-------------------------------------------------------------------------------------------------------
def one_screening_parameter_plot(
    x: np.ndarray[float],
    y: np.ndarray[float],
    labels: Iterable[str],
    unit_of_measure: str,
    scale: Literal['linear', 'discrete', 'log'],
    ax_lims: Tuple[Tuple[float], Tuple[float]],
    ax: matplotlib.axes._axes.Axes
) -> None:
    '''
    Add one subplot to the dashboard in the case the screening parameter is just one.

    x: (np.ndarray[float])
        The screened parameter to put on the x-axis in this subplot.
        
    y: (np.ndarray[float])
        The feature to put on the x-axis in this subplot.
        It can be either a 1D array, or a 2D array (e.g., mean + std deviation).

    labels: (Iterable[str])
        The names of the features on x and y axes.
    
    unit_of_measure: (str)
        The unit of measure associated to the x-axis.
    
    scale: (Literal['linear', 'discrete', 'log'])
        A string that specifies the scale to use on x-axis.

    ax_lims: (Tuple[Tuple[float], Tuple[float]])
        The x-axis limit for the current plot.

    ax: (matplotlib.axes._axes.Axes)
        The ax object associated to the subplot to add.
    '''

    # Get parameter to put on x-axis
    x_name = labels[0]
    x_unit = unit_of_measure

    # Get indicator to put on y-axis
    y_name = labels[1]
    if y_name == 'delta_sphericity':
        y_name = r"$\Delta$_sphericity"
    elif y_name == 'IoU':
        y_name = "mean(1 - IoU)"
    elif y_name == 'IoU_derivative':
        y_name = "mean(1 - IoU)_derivative"
    elif y_name == 'volume_loss':
        y_name = "% volume_loss"
    
    if len(y.shape) > 1:
        ax.scatter(x, y[:, 0], s=100)
        ax.errorbar(x, y[:, 0], yerr=y[:, 1], capsize=5, fmt='--', ecolor='grey', markersize=2)
    else:
        ax.scatter(x, y, s=100)
        ax.plot(x, y, linestyle='dashed')

    # Set labels with unit of measure
    if x_unit:
        xlab = f"{x_name.replace('_', ' ').title()} ({x_unit})"
    else:
        xlab = f"{x_name.replace('_', ' ').title()}"
        
    # Set x-axis stuff
    ax.set_xlabel(xlab, fontsize=18)
    if scale == 'log':
        ax.set_xscale('log')
    elif scale == 'discrete':
        discrete_x_values = np.unique(x)
        ax.xaxis.set_major_locator(matplotlib.ticker.FixedLocator(discrete_x_values))
        ax.xaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        ax.xaxis.set_minor_formatter(ScalarFormatter(useMathText=True))
    
    # Set y-axis stuff
    ylab = f"{y_name.replace('_', ' ').title()}"
    ax.set_ylabel(ylab, fontsize=18)

    # Set axes limits
    ax.set_xlim(left=ax_lims[0][0], right=ax_lims[0][1])
    ax.set_ylim(bottom=ax_lims[1][0], top=ax_lims[1][1])
#-------------------------------------------------------------------------------------------------------



#-------------------------------------------------------------------------------------------------------
def one_param_failed_plot(
    failed: np.ndarray[float],
    outliers: np.ndarray[float],
    interpenetration: np.ndarray[float],  
    label: str,
    unit_of_measure: str,
    scale: Literal['linear', 'discrete', 'log'],
    ax: matplotlib.axes._axes.Axes
) -> Tuple[Tuple[float], Tuple[float]]:
    """
    Add a subplot that report failed, interpenetrating and outlying simulation runs 
    in the case of one single screening parameter.

    Parameters:
    -----------

    failed: (np.ndarray[float])
        The array of screened parameters correspondent to failed simulation runs.

    outliers: (np.ndarray[float])
        The array of screened parameters correspondent to outlying simulation runs.

    interpenetration: (np.ndarray[float])
        The array of screened parameters correspondent to simulation runs with interpenetration.

    label: (str)
        The name of the parameter to show on x-axis.
    
    unit_of_measure: (str)
        The unit of measure associated to the x-axis.
    
    scale: (Literal['linear', 'discrete', 'log'])
        A string that specifies the scale to use on x-axis.

    ax: (matplotlib.axes._axes.Axes)
        The ax object associated to the subplot to add.

    Returns:
    --------
    
    (Tuple[Tuple[float], Tuple[float]]):
        The axes limits for the current plot.
    """

    # Get mid point on y-axis
    y_min, y_mid, y_max = 0, 0.5, 1

    # Plot points of failed runs as 'x'
    ax.scatter(
        failed, [y_mid]*len(failed), s=100, 
        marker='x', color='red', label='Failed'
    )

    # Plot outlying points as 'o'
    ax.scatter(
        outliers, [y_mid]*len(outliers), 
        facecolors='none', edgecolors='green', 
        linewidths=1.5, label='Outlier', s=100
    )

    # Plot interpenetrating points as a triangle
    ax.scatter(
        interpenetration, [y_mid]*len(interpenetration), 
        marker='^', facecolors='none', edgecolors='orange', 
        linewidths=1.5, label='Interpenetrating', s=100
    )

    # Add legend for failed, outlying and interpenetrating points
    ax.legend(loc="best")

    # Set y-axis limits
    ax.set_ylim([y_min, y_max])

    # Add unit of measures to axes labels
    if unit_of_measure:
        xlab = f"{label.replace('_', ' ').title()} ({unit_of_measure})"
    else:
        xlab = f"{label.replace('_', ' ').title()}"
    
    # Set axes labels
    ax.set_xlabel(xlab, fontsize=18)

    # Set x-axis stuff
    if scale == 'log':
        ax.set_xscale('log')
    elif scale == 'discrete':
        discrete_x_values = np.unique(np.concatenate([failed, interpenetration, outliers]))
        ax.xaxis.set_major_locator(matplotlib.ticker.FixedLocator(discrete_x_values))
        ax.xaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        ax.xaxis.set_minor_formatter(ScalarFormatter(useMathText=True))

    # Add horizontal line for y_mid
    ax.axhline(y=y_mid, color='gray', linestyle='dashed', linewidth=1)

    return ax.get_xlim(), ax.get_ylim()
#-------------------------------------------------------------------------------------------------------



#-------------------------------------------------------------------------------------------------------
def two_screening_parameters_plot(
    x: np.ndarray[float],
    y: np.ndarray[float],
    z: np.ndarray[float],
    labels: Iterable[str],
    units_of_measure: Iterable[str],
    scales: Iterable[Literal['linear', 'discrete', 'log']],
    ax_lims: Tuple[Tuple[float], Tuple[float]],
    ax: matplotlib.axes._axes.Axes
) -> None:
    '''
    Add one subplot to the dashboard in the case of exactly 2 screening parameters.

    x: (np.ndarray[float])
        The screened parameter to put on the x-axis in this subplot.
        
    y: (np.ndarray[float])
        The feature/parameter to put on the y-axis in this subplot.
        It can be either a 1D array, or a 2D array (e.g., mean + std deviation).
    
    z: (np.ndarray[float])
        If there are at least 2 screening parameters, this is the array of values of a feature
        for every parameters combination.

    labels: (Iterable[str])
        The names of the features associated to x, y, z arrays.
    
    units_of_measure: (Iterable[str])
        The units of measure associated to the axes.
    
    scales: (Iterable[Literal['linear', 'discrete', 'log']])
        A pair of strings that specify the scale to use on each axis.
    
    ax_lims: (Tuple[Tuple[float], Tuple[float]])
        The axes limits for the current plot.

    ax: (matplotlib.axes._axes.Axes)
        The ax object associated to the subplot to add.
    '''

    # Plot data
    if len(z.shape) > 1:
        sc = ax.scatter(x, y, c=z[:, 0], cmap='viridis', s=100)
        # z_std_norm = (z[:, 1] - min(z[:, 1])) / (max(z[:, 1]) -  min(z[:, 1])) * 500
        # ax.scatter(x, y, s=z_std_norm, alpha=0.7, facecolors='none', edgecolors='grey', linewidths=1.5, label='Std Error')
    else:
        sc = ax.scatter(x, y, c=z, cmap='viridis', s=100)

    # Get names for axes labels
    x_name = labels[0]
    x_unit = units_of_measure[0]
    y_name = labels[1]
    y_unit = units_of_measure[1]
    z_name = labels[2]
    if z_name == 'delta_sphericity':
        z_name = r"$\Delta$_sphericity"
    elif z_name == 'IoU':
        z_name = "mean(1 - IoU)"
    elif z_name == 'IoU_derivative':
        z_name = "mean(1 - IoU)_derivative"
    elif z_name == 'volume_loss':
        z_name = "% volume_loss"

    # Add color bar
    cbar = plt.colorbar(sc, ax=ax)
    cbar.set_label(z_name.replace('_', ' ').title(), fontsize=14) 

    # Add unit of measures to axes labels
    if x_unit:
        xlab = f"{x_name.replace('_', ' ').title()} ({x_unit})"
    else:
        xlab = f"{x_name.replace('_', ' ').title()}"
    
    if y_unit:
        ylab = f"{y_name.replace('_', ' ').title()} ({y_unit})"
    else:
        ylab = f"{y_name.replace('_', ' ').title()}"

    # Set axes labels
    ax.set_xlabel(xlab, fontsize=10)
    ax.set_ylabel(ylab, fontsize=10)

    # Set x-axis stuff
    if scales[0] == 'log':
        ax.set_xscale('log')
    elif scales[0] == 'discrete':
        discrete_x_values = np.unique(x)
        ax.xaxis.set_major_locator(matplotlib.ticker.FixedLocator(discrete_x_values))
        ax.xaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        ax.xaxis.set_minor_formatter(ScalarFormatter(useMathText=True))
    
    # Set y-axis stuff
    if scales[1] == 'log':
        ax.set_yscale('log')
    elif scales[1] == 'discrete':
        discrete_y_values = np.unique(y)
        ax.yaxis.set_major_locator(matplotlib.ticker.FixedLocator(discrete_y_values))
        ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        ax.yaxis.set_minor_formatter(ScalarFormatter(useMathText=True))

    # Set axes limits
    ax.set_xlim(left=ax_lims[0][0], right=ax_lims[0][1])
    ax.set_ylim(bottom=ax_lims[1][0], top=ax_lims[1][1])

    # Add horizontal lines for unique 'y' values
    unique_y = np.unique(y)
    for y_value in unique_y:
        ax.axhline(y=y_value, color='gray', linestyle='dashed', linewidth=1)
#-------------------------------------------------------------------------------------------------------



#-------------------------------------------------------------------------------------------------------
def two_params_failed_plot(
    failed: np.ndarray[float],
    outliers: np.ndarray[float],
    interpenetration: np.ndarray[float],  
    labels: Iterable[str],
    units_of_measure: Iterable[str],
    scales: Iterable[Literal['linear', 'discrete', 'log']],
    ax: matplotlib.axes._axes.Axes
) -> Tuple[Tuple[float], Tuple[float]]:
    """
    Add a subplot that report failed, interpenetrating and outlying simulation runs 
    in the case of more than 2 screening parameters.

    Parameters:
    -----------

    failed: (np.ndarray[float])
        The (x,y) pairs of the screened parameters correspondent to failed simulation runs.

    outliers: (np.ndarray[float])
        The (x,y) pairs of the screened parameters correspondent to outlying simulation runs.

    interpenetration: (np.ndarray[float])
        The (x,y) pairs of the screened parameters correspondent to simulation runs with interpenetration.

    labels: (Iterable[str])
        The names of the features associated to x, y arrays.
    
    units_of_measure: (Iterable[str])
        The units of measure associated to the axes.
    
    scales: (Iterable[Literal['linear', 'discrete', 'log']])
        A pair of strings that specify the scale to use on each axis.

    ax: (matplotlib.axes._axes.Axes)
        The ax object associated to the subplot to add.

    Returns:
    --------

    (Tuple[Tuple[float], Tuple[float]]):
        The axes limits for the current plot.
    """

    data = pd.DataFrame({
        'x': np.concatenate([failed[:, 0], outliers[:, 0], interpenetration[:, 0]], axis=0),
        'y': np.concatenate([failed[:, 1], outliers[:, 1], interpenetration[:, 1]], axis=0),
        'subset': ["fail"] * len(failed) + ["out"] * len(outliers) + ["inter"] * len(interpenetration)
    })

    # Plot points of failed runs as 'x'
    if len(failed) > 0:
        subset_data = data[data["subset"] == "fail"]
        ax.scatter(subset_data["x"], subset_data["y"], s=100, marker='x', color='red', label='Failed')

    # Plot outlying points as 'o'
    if len(outliers) > 0:
        subset_data = data[data["subset"] == "out"]
        ax.scatter(subset_data["x"], subset_data["y"], s=100, facecolors='none', 
                   edgecolors='green', linewidths=1.5, label='Outlier')

    # Plot interpenetrating points as a triangle
    if len(interpenetration) > 0:
        subset_data = data[data["subset"] == "inter"]
        ax.scatter(subset_data["x"], subset_data["y"], s=100,  marker='^', 
                   facecolors='none', edgecolors='orange', linewidths=1.5, label='Interpenetrating')

    # Add legend for failed, outlying and interpenetrating points
    ax.legend(loc="best")

    # Add unit of measures to axes labels
    if units_of_measure[0]:
        xlab = f"{labels[0].replace('_', ' ').title()} ({units_of_measure[0]})"
    else:
        xlab = f"{labels[0].replace('_', ' ').title()}"
    
    if units_of_measure[1]:
        ylab = f"{labels[1].replace('_', ' ').title()} ({units_of_measure[1]})"
    else:
        ylab = f"{labels[1].replace('_', ' ').title()}"

    # Set axes labels
    ax.set_xlabel(xlab, fontsize=10)
    ax.set_ylabel(ylab, fontsize=10)

    # Set x-axis stuff
    if scales[0] == 'log':
        ax.set_xscale('log')
    elif scales[0] == 'discrete':
        discrete_x_values = np.unique(data["x"])
        ax.xaxis.set_major_locator(matplotlib.ticker.FixedLocator(discrete_x_values))
        ax.xaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        ax.xaxis.set_minor_formatter(ScalarFormatter(useMathText=True))
    
    # Set y-axis stuff
    if scales[1] == 'log':
        ax.set_yscale('log')
    elif scales[1] == 'discrete':
        discrete_y_values = np.unique(data["y"])
        ax.yaxis.set_major_locator(matplotlib.ticker.FixedLocator(discrete_y_values))
        ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        ax.yaxis.set_minor_formatter(ScalarFormatter(useMathText=True))

    # Add horizontal lines for unique 'y' values
    unique_y = data['y'].unique()
    for y_value in unique_y:
        ax.axhline(y=y_value, color='gray', linestyle='dashed', linewidth=1)

    return ax.get_xlim(), ax.get_ylim()
#-------------------------------------------------------------------------------------------------------



#-------------------------------------------------------------------------------------------------------
def many_screening_parameters_plot(
    x: np.ndarray[float],
    y: np.ndarray[float],
    z: np.ndarray[float],
    labels: Iterable[str],
    units_of_measure: Iterable[str],
    scales: Iterable[Literal['linear', 'discrete', 'log']],
    ax_lims: Tuple[Tuple[float], Tuple[float]],
    ax: matplotlib.axes._axes.Axes
) -> None:
    '''
    Add one subplot to the dashboard in the case of more than 2 screening parameters.

    x: (np.ndarray[float])
        The screened parameter to put on the x-axis in this subplot.
        
    y: (np.ndarray[float])
        The feature/parameter to put on the y-axis in this subplot.
        It can be either a 1D array, or a 2D array (e.g., mean + std deviation).
    
    z: (np.ndarray[float])
        If there are at least 2 screening parameters, this is the array of values of a feature
        for every parameters combination.

    labels: (Iterable[str])
        The names of the features associated to x, y, z arrays.
    
    units_of_measure: (Iterable[str])
        The units of measure associated to the axes.
    
    scales: (Iterable[Literal['linear', 'discrete', 'log']])
        A pair of strings that specify the scale to use on each axis.
    
    ax_lims: (Tuple[Tuple[float], Tuple[float]])
        The axes limits for the current plot.

    ax: (matplotlib.axes._axes.Axes)
        The ax object associated to the subplot to add.
    '''

    # Set axes scales (with swarmplot this must be done in advance)
    if scales[0] == 'log':
        ax.set_xscale('log')
    if scales[1] == 'log':
        ax.set_yscale('log')


    # Plot data
    palette1 = sns.color_palette("viridis", as_cmap=True)
    sns.swarmplot(x=x, y=y, hue=z, palette=palette1, size=3, ax=ax)

    # Get variable names 
    x_name = labels[0]
    x_unit = units_of_measure[0]
    y_name = labels[1]
    y_unit = units_of_measure[1]
    z_name = labels[2]
    if z_name == 'delta_sphericity':
        z_name = r"$\Delta$_sphericity"
    elif z_name == 'IoU':
        z_name = "mean(1 - IoU)"
    elif z_name == 'IoU_derivative':
        z_name = "mean(1 - IoU)_derivative"
    elif z_name == 'volume_loss':
        z_name = "% volume_loss"
    
    # Add color bar
    sm = plt.cm.ScalarMappable(cmap=palette1)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_label(z_name.replace('_', ' ').title())

    # Hide legend
    ax.legend([],[], frameon=False)

    # Add unit of measures to axes labels
    if x_unit:
        xlab = f"{x_name.replace('_', ' ').title()} ({x_unit})"
    else:
        xlab = f"{x_name.replace('_', ' ').title()}"
    
    if y_unit:
        ylab = f"{y_name.replace('_', ' ').title()} ({y_unit})"
    else:
        ylab = f"{y_name.replace('_', ' ').title()}"

    # Set axes labels
    ax.set_xlabel(xlab, fontsize=10)
    ax.set_ylabel(ylab, fontsize=10)

    # Set axes limits
    ax.set_xlim(left=ax_lims[0][0], right=ax_lims[0][1])
    ax.set_ylim(bottom=ax_lims[1][0], top=ax_lims[1][1])

    # Set x-axis stuff
    if scales[0] == 'log':
        ax.set_xticks([])
    elif scales[0] == 'discrete':
        discrete_x_values = np.unique(x)
        ax.xaxis.set_major_locator(matplotlib.ticker.FixedLocator(discrete_x_values))
        ax.xaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        ax.xaxis.set_minor_formatter(ScalarFormatter(useMathText=True))
    
    # Set y-axis stuff
    if scales[1] == 'log':
        ax.set_yticks([])
    elif scales[1] == 'discrete':
        discrete_y_values = np.unique(y)
        ax.yaxis.set_major_locator(matplotlib.ticker.FixedLocator(discrete_y_values))
        ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        ax.yaxis.set_minor_formatter(ScalarFormatter(useMathText=True))

    # Add horizontal lines for unique 'y' values
    unique_y = y.unique()
    for y_value in unique_y:
        ax.axhline(y=y_value, color='gray', linestyle='dashed', linewidth=1)
#-------------------------------------------------------------------------------------------------------



#-------------------------------------------------------------------------------------------------------
def many_params_failed_plot(
    failed: np.ndarray[float],
    outliers: np.ndarray[float],
    interpenetration: np.ndarray[float],  
    labels: Iterable[str],
    units_of_measure: Iterable[str],
    scales: Iterable[Literal['linear', 'discrete', 'log']],
    ax: matplotlib.axes._axes.Axes
) -> Tuple[Tuple[float], Tuple[float]]:
    """
    Add a subplot that report failed, interpenetrating and outlying simulation runs 
    in the case of more than 2 screening parameters.

    Parameters:
    -----------

    failed: (np.ndarray[float])
        The (x,y) pairs of the screened parameters correspondent to failed simulation runs.

    outliers: (np.ndarray[float])
        The (x,y) pairs of the screened parameters correspondent to outlying simulation runs.

    interpenetration: (np.ndarray[float])
        The (x,y) pairs of the screened parameters correspondent to simulation runs with interpenetration.

    labels: (Iterable[str])
        The names of the features associated to x, y arrays.
    
    units_of_measure: (Iterable[str])
        The units of measure associated to the axes.
    
    scales: (Iterable[Literal['linear', 'discrete', 'log']])
        A pair of strings that specify the scale to use on each axis.

    ax: (matplotlib.axes._axes.Axes)
        The ax object associated to the subplot to add.

    Returns:
    --------

    (Tuple[Tuple[float], Tuple[float]]):
        The axes limits for the current plot.
    """
    
    # Get dataset for plotting
    data = pd.DataFrame({
        'x': np.concatenate([failed[:, 0], outliers[:, 0], interpenetration[:, 0]]),
        'y': np.concatenate([failed[:, 1], outliers[:, 1], interpenetration[:, 1]]),
        'subset': ["fail"] * len(failed) + ["out"] * len(outliers) + ["inter"] * len(interpenetration)
    })

    # Set axes scales (with swarmplot this must be done in advance)
    if scales[0] == 'log':
        ax.set_xscale('log')
    if scales[1] == 'log':
        ax.set_yscale('log')        

    # Overlap failed, outlying and interpenetrating runs
    palette2 = {
        'fail': to_rgba('red', 1),
        'out': to_rgba('green', 1),
        'inter': to_rgba('orange', 1)
    }
    sns.swarmplot(x='x', y='y', data=data, hue='subset', palette=palette2, size=3, marker='D', ax=ax)

    # Add legend for failed, outlying and interpenetrating points
    handles, lbls = ax.get_legend_handles_labels()
    runs_legend = ax.legend(handles=handles[-3:], labels=lbls[-3:])
    ax.add_artist(runs_legend)

    # Add unit of measures to axes labels
    if units_of_measure[0]:
        xlab = f"{labels[0].replace('_', ' ').title()} ({units_of_measure[0]})"
    else:
        xlab = f"{labels[0].replace('_', ' ').title()}"
    
    if units_of_measure[1]:
        ylab = f"{labels[1].replace('_', ' ').title()} ({units_of_measure[1]})"
    else:
        ylab = f"{labels[1].replace('_', ' ').title()}"

    # Set axes labels
    ax.set_xlabel(xlab, fontsize=10)
    ax.set_ylabel(ylab, fontsize=10)

    # Set x-axis stuff
    if scales[0] == 'log':
        ax.set_xticks([])
    elif scales[0] == 'discrete':
        discrete_x_values = np.unique(data["x"])
        ax.xaxis.set_major_locator(matplotlib.ticker.FixedLocator(discrete_x_values))
        ax.xaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        ax.xaxis.set_minor_formatter(ScalarFormatter(useMathText=True))
    
    # Set y-axis stuff
    if scales[1] == 'log':
        ax.set_yticks([])
    elif scales[1] == 'discrete':
        discrete_y_values = np.unique(data["y"])
        ax.yaxis.set_major_locator(matplotlib.ticker.FixedLocator(discrete_y_values))
        ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        ax.yaxis.set_minor_formatter(ScalarFormatter(useMathText=True))

    # Add horizontal lines for unique 'y' values
    unique_y = data['y'].unique()
    for y_value in unique_y:
        ax.axhline(y=y_value, color='gray', linestyle='dashed', linewidth=1)

    return ax.get_xlim(), ax.get_ylim()
#-------------------------------------------------------------------------------------------------------



#-------------------------------------------------------------------------------------------------------
def plot_dashboard(
    df: pd.DataFrame,
    param_names: Iterable[str],
    units_of_meas: Iterable[str],
    param_scales: Iterable[str],
    save_dir: str,
    show: Optional[bool] = True
) -> None:
    """
    Preprocess the input dataframe, extract data from it, and finally call the 
    function to generate the subplots for the dashboard.

    Parameters:
    -----------
    df: (pd.DataFrame)
        The input dataframe.

    param_names: (Iterable[str])
        The names of the screened parameters.

    units_of_meas: (Iterable[str])
        A list of units of measure asociated to `param_names`.

    param_scales: (Iterable[str])
        A list of axes scales asociated to `param_names`.

    save_dir: (str)
        The path to the directory for saving the plot. If `None` the plot is not saved.

    show: (Optional[bool] = True)
        If `True` the plot is shown.
    """
    
    # Extract indicators
    indicators = [
        column 
        for column in df.columns 
        if column not in param_names and column not in ['sim_id', 'interpenetration']
    ]
    
    # Find type of plots -> 3 possibilities:
    # - Only one parameter --> scatter plots: param on x-axis, indicator on y-axis
    # - Two parameters --> scatter plots: param on x, y axes, indicator as color/size of dots
    # - More than two parameters --> swarmplots: param on x, y axes, indicator as color/size of dots
    if len(param_names) == 1:
        plot_type = '1'
    elif len(param_names) == 2:
        plot_type = '2'
    else:
        plot_type = '3'
    
    # Find number of plots:
    # Schema --> each row for one indicator, number of columns equal to the parameters combination
    parameter_pairs = list(set(combinations(param_names, 2)))
    n_rows = len(indicators)
    n_cols = len(parameter_pairs)

    # Split columns with multiple values and find failed simulation runs
    clean_df, failed_df, outliers_df, interpen_df = preprocess_data(df, param_names)

    if plot_type == '1':

        # Create the plot figure
        fig = plt.figure(
            figsize=(24, 16),
            constrained_layout=True
        )
        fig.suptitle("Parameter Screening Dashboard", fontsize=30)

        # Get parameter to put on x-axis
        x_name = param_names[0]
        x = clean_df[x_name].values
        x_unit = units_of_meas[0]
        x_scale = param_scales[0]

        # Add subplot with failed, interpenetrating and outlying
        ax = fig.add_subplot(2, 3, 1)
        lim = one_param_failed_plot(
            failed=failed_df[x_name],
            outliers=outliers_df[x_name],
            interpenetration=interpen_df[x_name],
            label=x_name,
            unit_of_measure=x_unit,
            scale=x_scale,
            ax=ax
        )

        subplot_id = 2
        for i, indicator in enumerate(indicators):
            # Get the current axis object
            ax = fig.add_subplot(2, 3, subplot_id)
            subplot_id += 1

            # Get indicator to put on y-axis
            y_name = indicator
            if indicator + '_mean' in clean_df.columns:
                y = np.column_stack([
                    clean_df[y_name + '_mean'].to_numpy(), 
                    clean_df[y_name + '_std'].to_numpy()
                ])
            else:
                y = clean_df[y_name].to_numpy()

            one_screening_parameter_plot(
                x=x, y=y,
                unit_of_measure=x_unit,
                labels=[x_name, y_name],
                scale=x_scale,
                ax_lim=lim,
                ax=ax
            )
        
    elif plot_type == '2':

        # Create the plot figure
        fig = plt.figure(
            figsize=(24, 16),
            constrained_layout=True
        )
        fig.suptitle("Parameter Screening Dashboard", fontsize=30)

        # Get parameters to put on x and y axes
        x_name, y_name = parameter_pairs[0]
        x, y = clean_df[x_name], clean_df[y_name]
        fail_data = np.column_stack([failed_df[x_name], failed_df[y_name]]) if len(failed_df) else None
        out_data = np.column_stack([outliers_df[x_name], outliers_df[y_name]]) if len(outliers_df) else None
        inter_data = np.column_stack([interpen_df[x_name], interpen_df[y_name]]) if len(interpen_df) else None
        x_unit, y_unit = units_of_meas[0], units_of_meas[1]
        x_scale, y_scale = param_scales[0], param_scales[1]

        # Add subplot with failed, interpenetrating and outlying
        ax = fig.add_subplot(2, 3, 1)
        
        lims = two_params_failed_plot(
            failed=fail_data,
            outliers=out_data,
            interpenetration=inter_data,
            labels=(x_name, y_name),
            units_of_measure=(x_unit, y_unit), 
            scales=(x_scale, y_scale), 
            ax=ax
        )

        # Add subplots for all the indicators
        subplot_id = 2
        for i, indicator in enumerate(indicators):
            # Get the current axis object
            ax = fig.add_subplot(2, 3, subplot_id)
            subplot_id += 1

            # Get indicator values
            z_name = indicator
            if indicator + '_mean' in clean_df.columns:
                z = np.column_stack([
                    clean_df[indicator + '_mean'].to_numpy(), 
                    clean_df[indicator + '_std'].to_numpy()
                ])
            else:
                z = clean_df[indicator].to_numpy()

            two_screening_parameters_plot(
                x=x, y=y, z=z, 
                labels=(x_name, y_name, z_name),
                units_of_measure=(x_unit, y_unit), 
                scales=(x_scale, y_scale), 
                ax_lims=lims,
                ax=ax
            )
    
    elif plot_type == '3':

        # Create the plot figure
        fig = plt.figure(
            figsize=(n_cols * 8, n_rows * 5),
            constrained_layout=True
        )
        fig.suptitle("Parameter Screening Dashboard", fontsize=30)
        subfigs = fig.subfigures(n_rows + 1, 1)

        for j, parameter_pair in enumerate(parameter_pairs):
            # Get parameters to put on x and y axes
            x_name, y_name = parameter_pair
            x, y = clean_df[x_name], clean_df[y_name]
            fail_data = np.column_stack([failed_df[x_name], failed_df[y_name]]) if len(failed_df) else None
            out_data = np.column_stack([outliers_df[x_name], outliers_df[y_name]]) if len(outliers_df) else None
            inter_data = np.column_stack([interpen_df[x_name], interpen_df[y_name]]) if len(interpen_df) else None
            x_unit, y_unit = units_of_meas[0], units_of_meas[1]
            x_scale, y_scale = param_scales[0], param_scales[1]
            
            # Add plot for failed, interpenetrating and outlying runs
            subfig = subfigs[0]
            subplot_id = j + 1
            ax = subfig.add_subplot(1, n_cols, subplot_id)

            lims = many_params_failed_plot(
                failed=fail_data,
                outliers=out_data,
                interpenetration=inter_data,  
                labels=[x_name, y_name],
                units_of_measure=[x_unit, y_unit],
                scales=[x_scale, y_scale],
                ax=ax
            )

            for i, indicator in enumerate(indicators):
                subfig = subfigs[i+1] 
                # subfig.suptitle(f"{indicator.replace('_', ' ').title()}", fontsize=24)
                ax = subfig.add_subplot(1, n_cols, subplot_id)

                # Get indicator values
                z_name = indicator
                if indicator + '_mean' in clean_df.columns:
                    z = clean_df[indicator + '_mean'].values
                else:
                    z = clean_df[indicator].values

                many_screening_parameters_plot(
                    x=x, y=y, z=z, 
                    labels=(x_name, y_name, z_name),
                    units_of_measure=(x_unit, y_unit), 
                    scales=(x_scale, y_scale), 
                    ax_lims=lims,
                    ax=ax
                )

    # Save the current plot
    if save_dir:
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        existing_names = [
            fname for fname in os.listdir(save_dir)
            if fname.split("_")[-2] == plot_type
        ]
        if not existing_names:
            save_name = f"screening_dashboard_type_{plot_type}_v1.jpg"
        else:
            last_version = sorted(
                list(map(lambda x: int(x.split("_")[-1].split(".")[0][1:]), existing_names))
            )[-1]
            save_name = f"screening_dashboard_type_{plot_type}_v{last_version + 1}.jpg"
        plt.savefig(os.path.join(save_dir, save_name))
        
    # Show the plot
    if show:
        plt.show()
    else:
        plt.close()

    return
#-------------------------------------------------------------------------------------------------------



if __name__ == '__main__':

    ##### USER SECTION #####
    ### 1. SPECIFY INDICATORS TO INCLUDE IN THE DASHBOARD ###
    '''AVAILABLE INDICATORS: 
    'delta_sphericity', 'interpenetration', 'volume_loss', 'num_iterations', 'IoU', 'IoU_derivative
    NOTE: 'interpenetration' must always be included!
    '''
    indicators = ['delta_sphericity', 'interpenetration', 'volume_loss', 'num_iterations', 'IoU', 'IoU_derivative']

    ### 2. SPECIFY SCREENING PARAMETERS & RELATIVE AXES SCALES  ###
    '''AVAILABLE PARAMETERS:
    "damping_coefficient"
    "time_step"
    "contact_cutoff_adhesion"
    "contact_cutoff_repulsion"
    "bulk_modulus"
    "adherence_strength"
    "repulsion_strength"
    "surface_tension"
    "bending_modulus"
    "ecm_bulk_modulus"
    "ecm_adherence_strength"
    "ecm_repulsion_strength"
    "ecm_surface_tension"
    "ecm_bending_modulus"
    '''
    parameters = ['param_1', 'param_2', 'param_3']
    
    '''AVAILABLE SCALES: ['linear', 'log', 'discrete']
    Note: prefer discrete over linear when the number of values for a screening parameter is low
    '''
    axes_scales = ['linear', 'log', 'discrete']

    ### 3. SPECIFY PATH TO SCREENING ROOT DIRECTORY & FILE NAME OF SCREENING TABLE (CSV)
    root_dir = "path/to/screening/results/directory"
    screening_csv_file_name = "screening_table.csv"
    ##############################################################################################


    ##############################################################################################
    ### NOT FOR USER ###

    assert 'interpenetration' in indicators, "'interpenetration' must always be included in the indicators!"

    # Make sure to be in the './screening_analysis' working directory
    os.chdir(os.path.abspath(os.path.dirname(__file__)))

    # Get unit of measures from file './misc/param2unitofmeas.json'
    with open("./misc/param2unitofmeas.json", "r") as file:
        uom_dict = json.load(file)
    units_of_measure = [uom_dict[param] for param in parameters]

    #path to the csv table containing the parameters and the number of iterations and if it has crashed
    screening_table_path = os.path.join(root_dir, "summary_folder", screening_csv_file_name)
    #path to files computed from `run_data_collection.py`
    data_path = os.path.join(root_dir, "dashboard_files")
    #path where to save the plot
    save_plot_path = os.path.join(root_dir, "dashboard_plots")

    # Load data from files
    original_df = merge_dataframes(
        path_to_screening_table=screening_table_path,
        path_to_data_dir=data_path,
        param_names=parameters,
        ind_names=indicators
    )

    # Plot (and process) the data
    plot_dashboard(
        df=original_df,
        param_names=parameters,
        units_of_meas=units_of_measure,
        param_scales=axes_scales,
        save_dir=save_plot_path,
        show=True
    )

    
