# Contrast-based GED cleaning for ambulatory EEG
> ⭐ Main GED implementation: `src/ged_eeg_analysis/contrast_cleaning.py`


## Table of Contents
- [Introduction](#introduction)
- [Running the Example](#running-the-example)
- [Results from the Figs Directory](#results-from-the-figs-directory)
- [Contributing](#contributing)
- [License](#license)

## Introduction
Motion artifacts in EEG recordings, whether from subject movement, electrode displacement, or environmental interference, introduce noise that can obscure meaningful brain activity. This project focuses on developing an unsupervised method for decomposing multichannel EEG to isolate large amplitude artifacts such as motion.

## Running the Example

To run the example code, follow these steps:

1. **Clone the repository:**

    ```bash
    git clone git@github.com:SaharSattari/ged_eeg_analysis.git
    ```

2. **Navigate to the project directory:**

    ```bash
    cd ged_eeg_analysis
    ```

3. **Install the required dependencies:**

    ```bash
    pip install -r requirements.txt
    ```

4. **Execute the example script:**

    ```bash
    python3 example/ged_analysis_simulated_data.py
    ```

    The script will run the analysis and display the results in your console.

## Results from the Figs Directory
Below are some key figures illustrating the example results:

- **Simulated Data:**

  ![Simulated data](figs/Artifactual.png)

- **Post-GED Analysis:**

  ![Post-GED simulated data](figs/PostGED.png)

- **GED Components:**

  ![GED Components](figs/GED_Components.png)

- **First Component Time Series:**
  ![FirstComponentTimeSeries](figs/FirstComponentTimeSeries.png)

## Reference
Sattari, S., Virji-Babul, N., & Wu, L. C. (2025). Contrast-based artifact removal enables microstate analysis in ambulatory EEG. IEEE Transactions on Biomedical Engineering. DOI: [10.1109/TBME.2025.3630112]([https://doi.org/10.1109/TBME.2025.3630112](https://ieeexplore.ieee.org/abstract/document/11231101/?casa_token=_kscB377Vz8AAAAA:o9dVfpvHmiyyyriKd6LZbHQS0Ramd4RkQXmj17MxlBQpJXjlbLiqp-Hj3se84Clevvf0wIGpJPs)
