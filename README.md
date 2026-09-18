# Electroencephalogram---EEG-# Quantitative Electroencephalogram (EEG) Analysis & Emotion Classification

This repository contains the digital signal processing (DSP) codebase and comprehensive research documentation for a Biomedical Engineering study conducted at Tel Aviv University. The project investigates the isolation of neuro-oscillatory bands, the automated detection of electrolyte salt bridges, and the continuous decoding of internal emotional states (Happy vs. Sad) using predictive machine learning.

## Clinical & Engineering Highlights

* **Advanced Signal Pre-processing:** Implemented custom 1-40 Hz FIR bandpass filters to attenuate high-frequency environmental noise and isolate continuous Delta (0.5-4 Hz), Theta (4-8 Hz), Alpha (8-13 Hz), and Beta (13-30 Hz) brainwaves.
* **Automated Artifact Detection:** Developed a normalized Electrical Distance (ED) algorithm to scan high-density 14-channel electrode montages for localized low-impedance electrolyte salt bridges.
* **Independent Component Analysis (ICA):** Utilized EEGLAB and the ICLabel extension to statistically isolate true cortical generators from non-neural ocular (99.7% probability) and cardiac (96.4% probability) artifacts.
* **Predictive Emotion Classification:** Extracted Power Spectral Density (PSD) features using a modified Welch periodogram and trained a linear Support Vector Machine (SVM) to decode pre-movement neural markers of imagined emotions.

## Repository Architecture

* `EEG - Paz and Sharon.pdf`
  The complete final report containing physiological background, mathematical methodologies, spatial heatmaps, temporal waveforms, and extended clinical discussions.
* `part11.m`
  MATLAB script responsible for primary signal cleaning. It applies the FIR bandpass filters, calculates standard deviations for isolated EEG rhythms across different ocular states, and generates frequency domain spectrums via Fast Fourier Transform (FFT).
* `part22.m`
  Script executing the Electrical Distance algorithm to identify adjacent channel short-circuits. It plots absolute signal differences and outputs a normalized ED matrix.
* `part33.m`
  The core script for emotion classification. It applies 4th-order IIR Butterworth filters, extracts spectral power features from selected Independent Components (IC2 and IC3), and evaluates a 10-Fold Cross-Validated linear SVM classifier with dual ROC curves.

## Key Experimental Results

* **Ocular State Modulation:** The isolated Alpha rhythm demonstrated a significant physiological amplitude increase during eye closure (12.21 µV) compared to relaxed wakefulness (3.68 µV).
* **Electrode Topography Integrity:** The automated ED algorithm successfully identified the P8-T8 channel configuration as the primary short-circuit risk with a minimal normalized ED of 0.1771. 
* **Machine Learning Performance:** Following operating threshold calibration, the linear SVM emotion classifier successfully achieved a specificity of 79.17% and an operational Area Under the Curve (AUC) of 0.69 when decoding "Happy" vs. "Sad" internal states.

## How to Run

1. Clone the repository and ensure MATLAB with the Signal Processing Toolbox and Statistics and Machine Learning Toolbox is installed.
2. Ensure the raw CSV/TXT data files are located in the same working directory as the `.m` scripts.
3. Run `part11.m` to observe baseline drift correction, FFT transformations, and standard deviation calculations across isolated wavebands.
4. Run `part22.m` to generate the overlaid channel difference plots and output the normalized Electrical Distance screening table.
5. Run `part33.m` to execute the Welch PSD feature extraction, train the SVM, and generate the threshold-optimized ROC performance curves.
