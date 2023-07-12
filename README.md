# Seismic spectrogram analysis

The code provides a set of Python functions for importing, visualizing, and analyzing seismic data. The functions allow the user to import a seismic volume from a SEGY file, calculate correlations between neighboring traces, obtain the average trace, and perform spectrogram analysis on a seismic trace. The code is open source and can be modified to fit specific needs.

The code can be broken down into 5 main blocks, each of which is defined by a specific function.

* **Importing the SEGY file:** `import_segy_volume`<br>
Imports a SEGY volume from a specified file into Python's environment.<br>


* **Plotting correlation maps:** `correlation_heatmap`<br>
Plots a correlation map using a heatmap.<br>


* **Calculating the average trace from the seismic volume:** `average_trace`<br>
Calculates the average trace from neighboring traces around a central trace.<br>


* **Plotting seismic traces:** `trace_plot`<br>
Plots two seismic traces, allowing the user to compare any arbitrary traces.<br>


* **Performing spectrogram analysis:** `spectrogram_analysis`<br>
Performs spectrogram analysis of a given trace signal.<br>


## Neighboring Traces Analysis: Correlation Maps
<center>
  <img src="/images/02_correlation_map.png" alt="correlation_map" />
</center>

## Trace Analysis: Main Trace vs Average Trace 
<center>
  <img src="/images/03_trace_plot.png" alt="trace_plot" />
</center>

## Frequency and Phase Analysis: Trace Spectrogram
<center>
  <img src="/images/01_seismic_spectrogram.png" alt="trace_spectrogram" />
</center>
