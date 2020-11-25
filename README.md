![](https://travis-ci.com/uwescience/whaletracks.svg?branch=master)

# whaletracks
Analysis of ocean bottom seismometer (OBS) data to detect and track whales.

## Installation For Beginners:

Install [github desktop](https://desktop.github.com/) (or you can use git in the terminal if you prefer)

Install [Anaconda](https://www.anaconda.com/products/individual) python distribution for your operating system

Clone whaletracks to your local machine using Github Desktop 

If using Mac or Linux:
-   Open bash terminal and navigate to whaletracks directory using "cd whaletracks"
-   Setup python package by running "python setup.py develop"
-   Download dependencies by running "pip install whaletracks"

If using Windows:
-   Open Anaconda Navigator
-   Launch "CMD.exe Prompt"
-   In CMD.exe terminal, use "cd" command to navigate to whaletracks folder
-   Setup python package by running "python setup.py develop"
-   Download dependencies by running "pip install whaletracks"

Now you should be able to open any python code in the whaletracks directory using the Spyder IDE from the Anaconda Navigator. Use Spyder to edit any code parameters. 

It is recommended that you run these codes using a terminal, NOT Spyder. Save any code edits you made from Spyder, then open either a CMD.exe (Windows) or bash (Mac, Linux) terminal. In the terminal, enter 'ipython'. This will open an instance of python in your terminal. From there, run any python codes by entering 'run code_name.py'

## Codes of interest in "whaletracks/whaletracks/detection" folder
- main_detection.py runs automated spectrogram cross-correlation detections of fin or blue whale calls on specified OBSs and times. Currently set up to run example fin whale detections. Modify parameters in this code to run it on instruments and times of your choosing.
- BlueCall_manual_picker.py creates spectrograms for user identification and selection of blue whale A and B calls. Currently set up to run example blue whale detections. Modify parameters in this code to run it on instruments and times of your choosing.

- detect_calls.py defines functions used by main_detection.py (casual users of this code will not need to edit these)
- detect_manual.py defines functions used by BlueCall_manual_picker.py (casual users of this code will not need to edit these)

## Resources
- [IRIS metadata aggragator](http://ds.iris.edu/mda/7D/FC03D/?starttime=2014-09-07T00:00:00&endtime=2015-10-02T23:59:59) use this to find information on stations you wish to run the detectors on. 
- [IRIS rest interface](https://service.iris.edu/irisws/fedcatalog/1/) (includes querying stations by geographic coordinates)
- [Query map by lat/lon](http://ds.iris.edu/gmap/#maxlat=50&maxlon=-124&minlat=38&minlon=-132&network=*&drawingmode=box&planet=earth)



