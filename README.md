# Orbit Related Calculation

This project calculates the line of sight between an Optical Ground Station (OGS) located in the Netherlands and the Sentinel-2B satellite. It uses the Skyfield library to perform orbital calculations and matplotlib to visualize the results.

## Setup

### Prerequisites

- Python 3.7 or higher
- Conda (recommended) or pip for package management

### Installation

1. Clone the repository:

    ```sh
    git clone https://github.com/yourusername/OrbitRelatedCalculation.git
    cd OrbitRelatedCalculation
    ```

2. Create and activate the conda environment:

    ```sh
    conda create --name orbitAnalyisis python=3.8
    conda activate orbitAnalyisis
    ```

3. Install the required packages:

    ```sh
    conda install numpy matplotlib skyfield
    ```

    Or if you are using pip:

    ```sh
    pip install numpy matplotlib skyfield
    ```

## Usage

1. Ensure the conda environment is activated:

    ```sh
    conda activate orbitAnalyisis
    ```

2. Run the script:

    ```sh
    python OgsOrbitcalc.py
    ```

## Code Overview

- [OgsOrbitcalc.py](http://_vscodecontentref_/0): Main script that calculates the line of sight between the OGS and the Sentinel-2B satellite. It also generates various plots to visualize the results.

## Output

The script generates several plots:

1. Histogram of Maximum Elevations
2. Histogram of Pass Durations
3. Histogram of Duration Above Elevation Threshold
4. Scatter Plot of Azimuth Change
5. Scatter Plot of Azimuth Rate
6. Scatter Plot of Maximum Elevation vs. Time
7. Scatter Plot of Duration Above Elevation Threshold vs. Pass Index (Seconds)

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

- [Skyfield](https://rhodesmill.org/skyfield/) for orbital calculations
- [Matplotlib](https://matplotlib.org/) for plotting