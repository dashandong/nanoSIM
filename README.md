# nanoSIM

**nanoSIM** is a MATLAB-based implementation of image reconstruction for Structured Illumination Microscopy (SIM). It provides tools for processing SIM raw data, performing deconvolution, and reconstructing high-resolution images using advanced algorithms.

## Features

- **Two-step LR Deconvolution**: Implements a modified version of the [two-step Lucy-Richardson deconvolution algorithm](https://doi.org/10.1038/srep37149) with energy conservation.
- **Wide-Field, Deconvolved, Optical-Sectioned and SIM Image Reconstruction**: Processes raw SIM data to produce wide-field, deconvolved, and high-resolution SIM images.
- **Batch Processing**: Supports automated batch processing of multiple time points.
- **Energy Normalization**: Ensures consistent energy across imaging modalities for fair comparisons.

## File Descriptions

- **`nanoSIM.m`**: Main script for SIM image reconstruction. Handles raw data input, OTF processing, and iterative reconstruction.
- **`mydeconvlucy.m`**: Custom implementation of the Lucy-Richardson deconvolution without non-negative constraint.
- **`gen2Dotf.m`**: Generates 2D OTFs based on optical system parameters.
- **`imwritestack32.m`**: Writes 32-bit image stacks to TIFF files.

## Requirements

- MATLAB R2019b or later.
- Image Processing Toolbox.
- A system capable of handling large image stacks and performing FFT operations efficiently.

## Usage

1. Clone or download the repository.
2. Open MATLAB and navigate to the project directory.
3. Run `nanoSIM.m` to start the SIM reconstruction process.
4. Follow the prompts to select raw SIM data and OTF files.
5. Adjust parameters as needed to fine-tune the reconstruction.