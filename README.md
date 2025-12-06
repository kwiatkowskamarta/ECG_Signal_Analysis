# Fetal Breathing Movement Detection via Doppler Signals

## Project Overview
This project implements a non-invasive method for detecting Fetal Breathing Movements (FBM) using ultrasonic Doppler signals. Fetal breathing movements - rhythmic contractions of the diaphragm - are a critical indicator of the fetus's central nervous system development and overall physiological health.

The algorithm processes raw Doppler signals to identify characteristic breathing patterns, enabling the calculation of statistical parameters such as mean frequency and standard deviation to assess the regularity and intensity of fetal activity.

This project was developed at the Institute of Metrology and Biomedical Engineering, Warsaw University of Technology.

## Technical Details
The analysis is performed on quadrature Doppler signals with the following specifications:
* Signal Type: Complex Doppler signal (Quadrature)
* Sampling Frequency: 2000 Hz
* Ultrasound Frequency: 2 MHz
* Methodology:
    * Bandpass filtration to eliminate noise.
    * Amplitude extraction of the Doppler signal.
    * Frequency analysis to isolate pseudo-respiratory patterns.

## Repository Structure
The project files are organized as follows:

| File Name | Description |
| :--- | :--- |
| `ECG Signal Analysis.mlx` | Main Live Script. Contains the source code, visualizations, and step-by-step analysis logic. |
| `ECG Signal Analysis.m` | Plain Code. An exported version of the logic for standard MATLAB execution without formatting. |
| `RESP2.mat` | Dataset. The input quadrature signal data (2 columns) used for testing the algorithm. |
| `Project_article.pdf`| Project Article. The original scientific draft describing the methodology and results (in Polish). |
| `ECG Signal Analysis.pdf` | Full Report. A PDF export of the Live Script results and plots. |

## How to Run
1.  Clone the repository:
    ```bash
    git clone (https://github.com/kwiatkowskamarta/ECG_Signal_Analysis.git)
    ```
2.  Open MATLAB.
3.  Load the project: Open `ECG Signal Analysis.mlx`.
4.  Run the analysis: Ensure `RESP2.mat` is in the same directory. Run the script to generate the frequency analysis plots and statistics.

## Author
**Marta Kwiatkowska** - *Warsaw University of Technology*
