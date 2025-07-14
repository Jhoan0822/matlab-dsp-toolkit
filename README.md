# MATLAB Digital Signal Processing (DSP) Toolkit

![MATLAB](https://img.shields.io/badge/Made%20with-MATLAB-E97420?style=for-the-badge&logo=mathworks)
![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg?style=for-the-badge)

A comprehensive MATLAB application for the acquisition, processing, analysis, and visualization of digital signals. Developed as a project for the Digital Signals Processing course at Universidad Del Valle.

**Contributors:**
*   Alejandro Muñoz Gutierrez
*   Jhoan Eduardo Saldarriaga Serna
*   Luis Alfonso Manzano Bolaños
*   Junior Alejandro Gómez Alarcón

---

### ► Overview

This project provides a complete, user-friendly system for performing a wide variety of signal processing tasks. Built with **MATLAB App Designer**, the application features an intuitive graphical interface that allows users to import signals from multiple sources, apply mathematical operations and filters, and visualize the results in both the time and frequency domains. The system is designed with robust data validation and preprocessing to ensure the integrity of all operations.

![MAINGUI](https://github.com/user-attachments/assets/608908be-44d1-4746-9320-4f5598a662c4)


### ► Key Features

####  signal Inputs & Generation
*   **Audio Recording:** Record audio directly from a microphone in mono or stereo.
*   **File Import:** Load signals from `.wav` and text files.
*   **Hardware DAQ:** Acquire real-time data from an Arduino via a serial port.
*   **Synthetic Signal Generator:** Create a variety of standard signals, including:
    *   Sine, Cosine, Sinc
    *   Chirp (Swept-frequency)
    *   Sawtooth, Triangular, Square (Step)
    *   Ramp

![SI2](https://github.com/user-attachments/assets/0e08de81-f583-4271-a2e5-5c3c7552dd25)


####  signal Preprocessing & Validation
*   **Amplitude Normalization:** Normalize signals to ranges `[0, 1]`, `[-1, 1]`, or using standard (Z-score) normalization.
*   **Data Integrity:** Automatic detection of `NaN` or `Inf` values with an option for linear interpolation to correct the data.
*   **Mono/Stereo Conversion:** Seamlessly convert signals between mono and stereo formats, or merge two mono signals into a single stereo output.
*   **Signal Compatibility:** Automatically resamples and pads signals to match sampling rates and lengths before performing operations between them.
![SI1](https://github.com/user-attachments/assets/e16d9c96-ab17-41b5-b462-2842b90572a9)


#### ⏱️ Time-Domain Operations
*   **Basic Arithmetic:** Addition, Subtraction, Multiplication, Division, and Exponentiation between two signals.
*   **Convolution:** Perform convolution between two signals.
*   **Transformations:** Apply time shifting and time reflection (reversal) to a signal.
*   **Sampling Rate Modification:** Up-sample or down-sample signals by an integer factor.
![OPSIGNALS](https://github.com/user-attachments/assets/17281595-53df-44a3-9613-f4943ab0c869)
![OPERATIONS](https://github.com/user-attachments/assets/5de1daa1-2b2f-4c96-a1db-909cf0898917)


#### 📊 Filtering & Frequency-Domain Processing
*   **FIR & IIR Filtering:** Apply pre-defined Finite Impulse Response (FIR) and Infinite Impulse Response (IIR) filters. The application uses block-based convolution for FIR filters and difference equations for IIR filters.
*   **Fourier Analysis:** Visualize the signal's frequency content with:
    *   **Magnitude Spectrum:** View the amplitude of frequency components.
    *   **Phase Spectrum:** View the phase of frequency components.
*   **Spectrogram:** Analyze how the frequency content of a signal changes over time using a spectrogram view.
![FREQ](https://github.com/user-attachments/assets/4921d65e-2a87-45ef-8e3b-62c214b1b7fb)

![SPECTROGRAMA](https://github.com/user-attachments/assets/a04b6277-52a7-4f3c-8c13-2d0215bd3580)

![GRAPHS](https://github.com/user-attachments/assets/ed8854cb-42b8-4bec-af5d-9ea4bfcb7a4c)

#### 🎧 Audio Playback & Output
*   **Interactive Playback:** Play input or output audio signals.
*   **Playback Controls:** Includes stop, reverse playback, and an echo effect.
*   **Volume Control:** Adjust volume with constant, increasing, or decreasing envelopes.
*   **Save to File:** Save any processed signal as a `.wav` file with configurable bit depth.
  
![AUDIO](https://github.com/user-attachments/assets/72660c6f-6f30-4f06-82ef-eab0f1948eec)

### ► Technical Details

*   **Development Environment:** Built entirely in **MATLAB** using **App Designer** for the graphical user interface.
*   **Core Algorithms:** Implements fundamental DSP algorithms from scratch, including convolution by blocks and difference equation solvers, alongside leveraging MATLAB's powerful built-in functions for FFT, resampling, and audio I/O.
*   **Hardware Integration:** Communicates with an Arduino board using the `serialport` interface for data acquisition.

### ► How to Run the Application

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/your-username/matlab-dsp-toolkit.git
    ```
2.  **Open MATLAB:** Navigate to the cloned repository's directory in MATLAB.
3.  **Run the app:** Open the `Project_of_Dsp_App.mlapp` file. The App Designer window will open. Click the **Run** button in the App Designer toolbar.

### ► License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.
