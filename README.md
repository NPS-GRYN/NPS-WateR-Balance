# NPS WateR Balance Model
This repository houses the R version of NPS Water Balance model and the coupled NPS Water Balance - IHACRES model of streamflow. A comprehensive user guide can be found `[here](https://drive.google.com/file/d/18XgJSvIbnpL7oj97bJhvrwmZNjM4ErDY/view?usp=sharing)`. Below is a brief description of each folder in this repository. The primary scripts that run the model are also described (sub-scripts and supporting scripts are not included in this description).  

### Code
The Code folder contains supporting scripts that define functions used by both the water balance and IHACRES streamflow models, including functions to pull and process data, run the model, and calibrate the model.

### Data
The Data folder contains two files with data necessary for the code to run; when the code is run, additional intermediate data (such as climate projections or USGS stream gage data) is stored in the Data folder as CSV files.

### Streamflow
The Streamflow folder includes scripts to run a model of streamflow, either using user-provided parameters or calibration to a USGS stream gage, for a user-specified geographic region. 
<br><br>`Run_Streamflow_Model.R`: This is the primary script in the Streamflow folder. it has four parts: defining model parameters and setting initial values of variables, pulling historical USGS streamflow observations and meteorological forcing data, running the model (including optimization, if desired), and analyzing model results for historical and future time periods.  
<br>`Streamflow_Historical_Analysis.R`: The Streamflow_Historical_Analysis.R script provides visualizations and preliminary analyses of historical measured streamflow trends based on USGS gage data. The script can be run independently or in conjunction with the Run_Streamflow_Model.R script. 

### Validation
The Validation folder includes scripts and supporting files to validate whether the performance of the R model matches other versions of the NPS water balance model that have been developed, including an Excel version and a Python version. 
<br><br>`FrogRock_crosscheck.R`: This script allows users to run the water balance model and confirm that its output matches the original version of the model in Excel and a more recently developed version of the model in Python, using the test site of Frog Rock in Yellowstone National Park. 

### Water Balance
The Water Balance folder includes scripts to run a water balance model using user-provided parameters, with an optional step of calibrating the model to OpenET estimates of AET.
`Run_WB_Model.R`: This is the primary script in the Water Balance folder. It has four parts: defining model parameters and setting initial values of variables, pulling OpenET observations and meteorological forcing data, running the model (including optimization, if desired), and analyzing model results for historical and future time periods. Each of these parts are carried out by sub-scripts called in the main script. For more details, please refer to the user guide.

