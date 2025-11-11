# LéXPLORE PP Mooring

## Project Information

The data is collected within the frame of the [LeXPLORE project](https://wp.unil.ch/lexplore/) on Lake Geneva. 
The data is used and displayed on the [Datalakes website](https://www.datalakes-eawag.ch/).

## Sensors
The mooring line include different sensors that collect high-frequency data at various depth. 
It measures dissolved oxygen concentration (DO), Photosynthetic Active Radiation (PAR) and water temperature.

## Installation

- Clone the repository to your local machine using the command: 

 `git clone https://renkulab.io/gitlab/lexplore/ppmooring.git`
 
 Note that the repository will be copied to your current working directory.

- Use Python 3 and install the requirements with:

 `pip install -r requirements.txt`

 The python version can be checked by running the command `python --version`. In case python is not installed or only an older version of it, it is recommend to install python through the anaconda distribution which can be downloaded [here](https://www.anaconda.com/products/individual). 

## Usage

### Process new data

In order to process new data locally on your machine the file path needs to be adapted to your local file system. The following steps are therefore necessary: 

- Edit the `scripts/input_batch.bat` file. Change all the directory paths to match your local file system. This file contains all the file paths necessary to launch the batch scripts `runfile.bat`.

- Edit the `scripts/input_python.py` file. Change all the directory paths to match your local file system. This file contains all the directories where the python script outputs data to.

To process new data, place the data in the input directory which you specified in the `scripts/input_batch.bat` file. 
Double-clicking on the `runfile.bat` file will automaticall process all the data in the input directory and store the output 
in the directories specified in the `scripts/input_python.py` file. Note that the input files need to have the same 
directory structure as the already exisiting files in the `data/Level 0` folder. 

### Adapt/Extend data processing piepeline

The python script `scripts/main_ppmooring.py` defines the different processing steps while the python script `scripts/ppmooring.py` contains the python class thetis with all the corresponding 
class methods to process the data. To add a new processing or visualization step, a new class method can be created in the `ppmooring.py` file and the step can be added in `main_ppmooring.py` file.
Both above mentioned python scripts are independent of the local file system.

## Data

The data can be found in the folder `data`. The data is structured as follows:

### Data Structure

- **Level 0**: Raw data collected from the different sensors.

- **Level 1A**: Raw data stored to NetCDF and csv file where attributes (such as sensors used, units, description of data, etc.) are added to the data.

- **Level 1B**: Column with quality flags are added to the Level 1A data. Quality flag "1" indicates that the data point didn't pass the quality checks and further investigation is needed, quality flag "0" indicates that no further investiagion is needed.

- **Level 2**: For each type of measurements the data from all the different depth levels is transformed into a temporal grid with a temporal spacing of 10 minutes and stored to a 
NetCDF file spanning over a time period of 10 days.

## Quality assurance

Two levels of quality checks are included in ppmooring data. The first level of quality checks includes range validation, data type checking and flagging missing data. The second level combines sucessively various time-series smoothing and outlier detection methods.

#### Water temperature
**Quality check level 1:**
- Remove data during maintenance periods (cross validation between temperature and PAR sensors)
- Data range validation
- Missing data

**Quality check level 2:**
Successive outlier detections methods 

- kmeans clustering
- convolutional smoothing with various window type
- quantile ranges
- relevance variation rates of water temperature

#### Dissolved oxygen
**Quality check level 1:**
- Remove data during maintenance periods (cross validation between temperature and PAR sensors)
- Data range based on sensors validation
- Missing data

**Quality check level 2:**
Successive outlier detections methods 

- data range validation according to seasonality and depth using SHL2 long-term monitoring data as reference
- relevance variation rates of DO
- time-series smoothing based on Savitzky–Golay filter
- k means clustering
- convolutional smoothing with various window type

#### PAR
**Quality check level 1:**
- Remove data during maintenance periods (cross validation between temperature and PAR sensors)
- Data range validation
- Missing data

**Quality check level 2:**
Successive outlier detections methods 

- kmeans clustering
- relevance variation rates of PAR
- convolutional smoothing with various window type
