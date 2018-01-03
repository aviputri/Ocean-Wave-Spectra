# Ocean-Wave-Spectra
A tool to generate a spectrum of a set of sea water level (SWL) data, and then to create a synthetic SWL and classify the synthetic waves parameters.

There are Python and MATLAB scripts to run this project. The Python script is indeed much faster. 

MATLAB code:
- The main of the MATLAB scripts is processall.m, which will call the zerodown function, to analyze the individual waves out of the SWL data by zerodown-crossing.
- In the zerodown.m, we are using lagrpol function (lagrpol.m) to use Lagrange interpolation on the SWL data gaps, and then the individual waves will be sorted to get the wave parameters.
- sintetis.m is used to create the wave spectrum (calling the spectrum.m in the script), synthetic SWL (ema.m), and then to repeat again the zerodown.m to identify and sort the individual waves.

Python code:
- All the codes are within 1 script (main.py). In the main process, we call through the functions at the bottom of the file.
- Functions can be called separately, ie for just zerodown_data or synthetic, type 'import main' in python interpreter and enter the functions manually.
- deleteblank.py is used to delete lines of blank data (-) if there is any
- break.py is used to separate the data into hourly files
