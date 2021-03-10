### CH4cast_v1 user guide
#Step 1 - DOWNLOAD CH4cast_v1 code via terminal in Rstudio

1: Open a fresh R Studio window without any existing scripts

2: Locate the terminal tab (next to console tab on bottom left) and select it

3: In terminal, make a FOLDER anywhere on your local PC that is called --> CH4cast_v1

      For example:
      mkdir /Users/Owner/Desktop/CH4cast_v1
      
4: In terminal, now make the new Folder your directory using the command "cd ".

      For example:
      cd /Users/Owner/Desktop/CH4cast_v1
      
5: In terminal, run the following command to download the CH4cast code from Github. 

      For example:
      git clone https://github.com/ryanmclake/CH4cast_v1.git
      
6: After the cloning finishes, there should be a new folder in your working directory named "CH4cast_v1"

Step - 2 Run three forecasts of CH4 ebullition from summer 2019

7: Navigate to the new folder called CH4cast_v1 under the "files" panel in Rstudio

8: Go to the scripts folder and OPEN ALL OF THE SCRIPTS starting at 01 and then to 04

9: In the 01_DataCompile script, find the "Source" key on the top right --> Click "Source" and then let the script run until it is finished

Note --> this may take a while as it is downloading published data from EDI. 

10: When the script is complete, select the 02_ForecastModelsJAGS_AllTraps

11: In the 02_ForecastModelsJAGS_AllTraps script, find the "Source" key on the top right --> Click "Source" and then let the script run until it is finished

Note --> this may take a while as it is running JAGS iteratively over multiple dates

12: When 02_ForecastModelsJAGS_AllTraps is complete, select the tab 03_PersistenceNullModel

13: In the 03_PersistenceNullModel script, find the "Source" key on the top right --> Click "Source" and then let the script run until it is finished

Note --> this may take a while as it is running JAGS iteratively over multiple dates

10: When 03_PersistenceNullModel is complete, select the tab 04_Analysis_Figures

9: In the 04_Analysis_Figures script, find the "Source" key on the top right --> Click "Source" and then let the script run until it is finished

Note --> this may take a while as it is writing a bunch of figures

10: Navigate to "forecast_output/figures" to find the plots from the three days that you just forecasted! 

### Congrats! You have just run three forecast cycles of CH4 ebullition rates in a reservoir
For questions or _construcitve_ feedback please email me at ryan333@vt.edu.
