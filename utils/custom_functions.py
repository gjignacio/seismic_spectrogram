# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 11:12:11 2023

@author: XJIG04
"""
import segyio
import numpy as np
from scipy.stats import pearsonr
from scipy import signal
import seaborn as sns
import matplotlib.pyplot as plt


class Trace:
    """
    A class to represent a seismic trace.

    Args:
        ilxl: A tuple of two integers representing the inline and crossline indices of the trace.
        seismic_volume_data: A 3-dimensional NumPy array containing the seismic volume data.
        seismic_volume_meta: A segyio.segy.spec object containing the seismic volume metadata.

    Raises:
        ValueError: If `ilxl` is not a tuple of two integers, or if `seismic_volume_data` is not a 
                3-dimensional NumPy array, or if `seismic_volume_meta` is not a segyio.segy.spec object.

    Attributes:
        ilxl: The inline and crossline indices of the trace.
        ilxl_idx: The index of the trace in the seismic volume data.
        data: The data for the trace.

    Methods:
        update_trace_idx: Updates the `ilxl_idx` attribute.
        update_trace_data: Updates the `data` attribute.
    """
    
    def __init__(self, ilxl, seismic_volume_data, seismic_volume_meta):
        if isinstance(ilxl, tuple) and len(ilxl) == 2 and all(isinstance(val, int) for val in ilxl):
            self.ilxl = ilxl
        else:
            raise ValueError("ilxl must be a tuple of two integers.")
        
        if isinstance(seismic_volume_data, np.ndarray) and seismic_volume_data.ndim == 3:
            self.seismic_volume_data = seismic_volume_data
        else:
            raise ValueError("seismic_volume_data must be a 3-dimensional NumPy array.")
        
        if isinstance(seismic_volume_meta, segyio.segy.spec):
            self.seismic_volume_meta = seismic_volume_meta
        else:
            raise ValueError("seismic_volume_meta must be a segyio.segy.spec object.")
        
        self.ilxl_idx = self.update_trace_idx()
        self.data = self.update_trace_data()
    
    def update_trace_idx(self):
        if self.ilxl is not None and self.seismic_volume_meta is not None:
            il_idx = np.where(self.seismic_volume_meta.ilines == self.ilxl[0])[0][0]
            xl_idx = np.where(self.seismic_volume_meta.xlines == self.ilxl[1])[0][0]
            return (il_idx, xl_idx)
        else:
            return None
        
    def update_trace_data(self):
        if self.ilxl is not None and self.seismic_volume_data is not None:
            trace_data = self.seismic_volume_data[self.ilxl_idx[0],self.ilxl_idx[1]]
            return trace_data
        else:
            return None

def import_segy_volume(file_path, strict=True):
    """
    Import a SEGY volume from the specified file.

    Args:
        file_path (str): Path to the SEGY file.
        strict (bool, optional): If True, enable strict mode for SEGY file reading.
                                 Defaults to True.

    Returns:
        tuple: A tuple containing the imported SEGY volume as a NumPy array and the metadata.

    Raises:
        FileNotFoundError: If the specified SEGY file does not exist.
        segyio.SEGYLibException: If an error occurs while reading the SEGY file.

    Notes:
        - The SEGY file should have the required headers and be in IBM float format or IEEE float format.
        - The imported volume will have shape (samples, inlines, crosslines) following the SEGY file format.
        - The metadata contains information about the geometry and other properties of the SEGY volume.

    Example:
        >>> file_path = "path/to/seismic.sgy"
        >>> volume, metadata = import_segy_volume(file_path)
        >>> print(volume.shape)  # Output: (samples, inlines, crosslines)
        >>> print(metadata)  # Output: {'inlines': [...], 'crosslines': [...], ...}
    """
  
    with segyio.open(file_path, "r", strict=strict) as segy:
        # volume = segy
        volume = segyio.tools.cube(segy)
        metadata = segyio.tools.metadata(segy)
        sample_rate = segyio.tools.dt(segy)
        
        
    print("SEGY file imported")
    print("IL: {}-{}".format(metadata.ilines.min(),metadata.ilines.max()))
    print("XL: {}-{}".format(metadata.xlines.min(),metadata.xlines.max()))
    return volume, metadata, sample_rate


def neighboring_traces_correlation(central_trace, seismic_volume_data, seismic_volume_meta, neighbors=5):
    """
    Calculate the correlation coefficients between a central trace and its neighboring traces.

    Parameters:
        central_trace (Trace): The central trace object.
        seismic_volume_data (ndarray): 3D array representing the seismic volume data.
        seismic_volume_meta (segyio.segy.spec): Metadata associated with the seismic volume data.
        neighbors (int, optional): Number of neighboring traces to consider. Default is 5.

    Returns:
        ndarray: 2D meshgrid array of inline (IL) values.
        ndarray: 2D meshgrid array of crossline (XL) values.
        ndarray: 2D array of correlation coefficients between the central trace and neighboring traces.
    """
    # Unpacking of the central trace inline/crossline values.
    il0, xl0 = central_trace.ilxl[0], central_trace.ilxl[1]
    
    # Generate inline and crossline arrays for neighboring traces
    ils = np.arange(il0 - neighbors, il0 + neighbors + 1, 1, dtype=int)
    xls = np.arange(xl0 - neighbors, xl0 + neighbors + 1, 1, dtype=int)
    IL, XL = np.meshgrid(ils, xls)
    CORR = np.zeros_like(IL, dtype=float)
    
    for il_idx, il in enumerate(ils):
        for xl_idx, xl in enumerate(xls):
            try:
                trace = Trace(ilxl=(int(il), int(xl)),
                              seismic_volume_data=seismic_volume_data,
                              seismic_volume_meta=seismic_volume_meta)
                
                corr_coef, _ = pearsonr(central_trace.data, trace.data)
                CORR[il_idx, xl_idx] = corr_coef
            except:
                pass
    
    correlation_heatmap(X=IL, Y=XL, Z=CORR, title='Pearson\'s Correlation')
    
    return IL, XL, CORR

def correlation_heatmap(X, Y, Z, title):
    """
    Plot a correlation map using a heatmap.

    Parameters:
        X (ndarray): 2D meshgrid representing the x-coordinates.
        Y (ndarray): 2D meshgrid representing the y-coordinates.
        Z (ndarray): 2D array representing the correlation values.
        title (str): Title of the correlation map.

    Returns:
        None
    """
    # Create a figure and axes
    fig, ax = plt.subplots(figsize=(12, 10))
    
    # Plot the heatmap
    sns.heatmap(Z, cmap='rocket_r', ax=ax, xticklabels=X[0], yticklabels=Y[:, 0], annot=True, fmt=".3f")
    
    # Set the axis labels
    ax.set_xlabel('INLINE')
    ax.set_ylabel('XLINE')
    
    # Set the title
    ax.set_title(title)
    
    # Show the plot
    plt.show()

def average_trace(central_trace, seismic_volume_data, seismic_volume_meta, neighbors=0):
    """
    Calculate the average trace from neighboring traces around a central trace.

    Parameters:
        central_trace (Trace): The central trace object.
        seismic_volume_data (ndarray): 3D array representing the seismic volume data.
        seismic_volume_meta (segyio.segy.spec): Metadata associated with the seismic volume data.
        neighbors (int, optional): Number of neighboring traces to consider. Default is 0.

    Returns:
        ndarray: The average trace computed from the neighboring traces.
    """
    # Unpacking of the central trace inline/crossline values.
    il0, xl0 = central_trace.ilxl[0], central_trace.ilxl[1]
    
    # Generate inline and crossline arrays for neighboring traces
    ils = np.arange(il0 - neighbors, il0 + neighbors + 1, 1, dtype=int)
    xls = np.arange(xl0 - neighbors, xl0 + neighbors + 1, 1, dtype=int)
    IL, XL = np.meshgrid(ils, xls)

    # Create the array to save the neighbor traces.
    traces_to_average = np.zeros((len(ils), len(xls), seismic_volume_data.shape[2]))
    
    for il_idx, il in enumerate(ils):
        for xl_idx, xl in enumerate(xls):
            try:
                trace = Trace(ilxl=(int(il), int(xl)),
                              seismic_volume_data=seismic_volume_data,
                              seismic_volume_meta=seismic_volume_meta)
                traces_to_average[il_idx, xl_idx, :] = trace.data
            except:
                pass
            
    average_trace = np.mean(traces_to_average, axis=(0, 1))
    print("- {} traces averaged.".format(len(ils) * len(xls)))
    
    return average_trace

def trace_plot(trace1=None, trace2=None, label1=None, label2=None, dt=0.002):
    """
    Plots two seismic traces.

    Parameters:
        trace1 (ndarray): The first trace to plot.
        trace2 (ndarray): The second trace to plot.
        label1 (str): The label for the first trace.
        label2 (str): The label for the second trace.
        dt (float): The sampling interval.

    Returns:
        None.

    """
    
    # Create a figure and axes
    fig, ax = plt.subplots(figsize=(12, 5))
    t=np.arange(0,2.3+dt,dt)
    
    #Plot the traces:
    ax.plot(t,trace1, label=label1, c="#070094", linewidth=1.5, zorder=0)
    if trace2 is not None:
        ax.plot(t,trace2, label=label2, c="red", linewidth=1, zorder=1)

    # Other customizations:
    ax.set_xlabel('Time [seconds]')
    ax.set_ylabel('Amplitude')
    plt.legend()

    # Show the plot
    plt.show()
    
def spectrogram_analysis(trace=None, dt=0.002, wl=101, ws=4):
    """
    Perform spectrogram analysis of a given trace signal.

    Parameters:
        trace (ndarray): 1D array representing the input trace signal.
        dt (float): Time interval between samples in seconds.
        wl (int): Window length for the spectrogram analysis.
        ws (int): Window step size for the spectrogram analysis.

    Returns:
        None
    """
    # Set up time variables for frequency and phase analysis
    fs = 1 / dt  # Sampling frequency [Hz]
    t = np.arange(0, len(trace) * dt, dt)

    # AMPLITUDE AND PHASE SPECTRUM CALCULATION

    # Parameters for the spectrogram
    noverlap = wl - ws  # Window overlapping

    # Calculate the amplitude spectrogram
    frequencies, times, spectrogram = signal.spectrogram(trace, fs=fs, window='hann', nperseg=wl, noverlap=noverlap)

    # Calculate the maximum value in the spectrogram
    max_value = np.max(spectrogram)

    # Convert the maximum value to dB
    max_dB = 10 * np.log10(max_value)

    spectrogram = 10 * np.log10(spectrogram)

    # Subtract the maximum value in dB from all values in the spectrogram
    spectrogram -= max_dB

    # Calculate the phase of the original signal
    phase_frequencies, phase_times, phase = signal.spectrogram(trace,
                                                              fs=fs,
                                                              window='hann',
                                                              nperseg=wl,
                                                              noverlap=noverlap,
                                                              mode="phase")

    # PLOT

    # Create subplots
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(10, 15), sharey=True)

    # Plot the time series in the first subplot
    ax1.plot(trace, t)
    ax1.set_xlim(-max(abs(trace)), max(abs(trace)))
    ax1.invert_yaxis()
    ax1.set_xlabel('Amplitude')
    ax1.set_ylabel('Time (s)')
    ax1.set_title('Signal')

    # Plot the amplitude spectrogram in the second subplot
    max_freq = 125
    freq_plot = ax2.pcolormesh(frequencies, times, spectrogram.T,
                               shading='auto',
                               cmap='gist_rainbow_r')
    freq_cbar = fig.colorbar(freq_plot, ax=ax2, orientation="vertical")
    freq_cbar.set_label('Amplitude (dB)')
    for hz in [25, 50, 75, 100]:
        ax2.axvline(x=hz, color='black', linestyle='--', linewidth=.5, alpha=.5)
    ax2.set_xlim(0, max_freq)
    ax2.set_xlabel('Frequencies (Hz)')
    ax2.set_title('Amplitude Spectrogram')

    # Plot the phase variation in the third subplot
    max_phase = 90
    phase_plot = ax3.contourf(phase_frequencies,
                              phase_times,
                              phase.T,
                              cmap='RdBu_r',
                              levels=np.arange(-max_phase, max_phase + 1, 1))

    for hz in [25, 50, 75, 100]:
        ax3.axvline(x=hz, color='red', linestyle='--', linewidth=.5, alpha=.25)

    phase_cbar = fig.colorbar(phase_plot, ax=ax3, orientation="vertical")
    cbar_ticks = [-max_phase, -max_phase/2, 0, max_phase/2, max_phase]
    cbar_ticks_labels = [str(tick) for tick in cbar_ticks]
    phase_cbar.set_ticks(cbar_ticks)
    phase_cbar.set_ticklabels(cbar_ticks_labels)
    phase_cbar.set_label('Phase (°)')
    ax3.set_xlim(0, max_freq)
    ax3.set_xlabel('Frequencies (Hz)')
    ax3.set_title('Phase Spectrogram')

    # Adjust margins to avoid overlapping
    fig.tight_layout()

    # Show the plot
    plt.show()

    # Save the plot as an image file
    fig.savefig("Output.png")
