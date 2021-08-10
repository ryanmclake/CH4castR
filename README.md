# CH4castR user guide

<a href="url"><img src = "images/NSF.png" align="top" height="200" width="200" ></a>
<a href="url"><img src = "images/FLARE.jpg" align="top" height="200" width="200" ></a>
<a href="url"><img src = "images/CH4cast.png" align="top" height="200" width="200" ></a>

-----


:busts_in_silhouette: Ryan McClure, Quinn Thomas, Mary Lofton, Whitney Woelmer, and Cayelan Carey

:busts_in_silhouette: Special thanks to: 

Questions?  :email: ryan333@vt.edu, rqthomas@vt.edu, cayelan@vt.edu

-----

## Motivation

Thank you for checking out NEON-forecast-code. Freshwater lakes globally are increasingly threatened as a result of rapidly changing land use and climate ([Carpenter et al., 2011](https://www.annualreviews.org/doi/abs/10.1146/annurev-environ-021810-094524)). In response, developing forecast workflows has has emerged as a powerful tool to predict future environmental conditions in lakes in order to make informed management decisions for safety, health, and conservation ([Carey et al., 2021](); [Baracchini et al., 2020](https://www.sciencedirect.com/science/article/pii/S0043135420300658); [Page et al., 2018](https://www.sciencedirect.com/science/article/pii/S0043135418300605)). However, the discipline of forecasting in lakes is still in the early stages of making forecasts that are robust and reproducible. As a result, there is a dire need for open-source forecast workflows that are broadly applicable to many different lake ecosystems and flexible to different datastreams and local needs.

Here, we applied the FLAREr forecasting system ([Thomas et al., 2020](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2019WR026138)) to six NEON lakes to test FLAREr's robustness and scalability to other sites. The NEON lakes serve as an exemplar case to test FLARE because they have reliable, open-source datastreams in which new data can be acquired at relatively low latencies (<1.5 months). The goal of our forecast scaling study was to show that FLAREr is scalable to other lake ecosystems and can produce robust forecasts of water temperatures up to 35-days into the future. Altogether, we hope this workflow is a first step to building a community of lake and reservoir forecast practitioners that develop reliable forecast workflows and make informed decisions for future lake conservation and management.


### Step 2 - Run forecasts of CH4 ebullition rates from summer 2019

7: Navigate to the new folder called CH4castR under the "files" panel in Rstudio

8: Go to the scripts folder and OPEN ALL OF THE SCRIPTS starting at 01 and then to 06

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
