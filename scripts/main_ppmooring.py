# -*- coding: utf-8 -*-
# Notes to improve: metadata from Level0 was not stored in Level1A. Not necessary for the moment.
# Import Packages
#Todo list : 13793/13793 [00:07<00:00, 1792.99it/s] and duplicated run...change...

import sys
import os
import yaml
import numpy as np
import fnmatch
import pandas as pd


# Local modification
root = "../renku_ppmooring/ppmooring"

#Delete .DS_Store
folderdel = os.path.join(root,"data")
import glob, os

for root1, dirs, files in os.walk(folderdel):
    for file in files:
        if file.endswith('.DS_Store'):
            path = os.path.join(root1, file)
            os.remove(path)


folder_script = os.path.join(root,"scripts")
sys.path.append(folder_script)

# Import ppmooring functions

# Important = change root in ppmooring class in ppmooring.py
from ppmooring import DO_data, PAR_data, Temperature_data

#### Main Script

#### Process level 0 - level1A.
# Input Directory
#input_dir = str(sys.argv[1]).replace('\\', '/')
input_dir = os.path.join(root,"data/Level0")
input_dirs = {
    "DO":input_dir+"/DO",
    "PAR":input_dir+"/RBR_PAR",
    "Temp":input_dir+"/Temperature"
}

# Read YAML file containing directories

with open(os.path.join(folder_script,"input_python.yaml"), "r") as f:
    directories = yaml.load(f, Loader=yaml.FullLoader)

# Make sure the directories exist
for directory in directories.values():
    if not os.path.exists(directory):
        os.makedirs(directory)

#1. QAQC PAR first #
# To determine the maintenance campaigns for all ppmooring sensors

## PAR ##
# Initalize object with its attributes
PAR_files = PAR_data()
# Get list of files to reaPAR_files.get_files_to_read_Level0(input_dirs["PAR"])

# Renew list of files to read
PAR_files.get_files_to_read_Level0(input_dirs["PAR"])

##Check and delete empty files
for infile in PAR_files.filestoread:
    fsize = os.stat(infile).st_size
    if fsize <= 1000: #delete files <1000 bytes
        os.remove(infile)

# Loop over all the files to read
for infile in PAR_files.filestoread:
    ## PPMooring PAR Data Level 0 --> Level 1A

    infile = infile.replace(os.sep,'/')

    # Create log files containing treated files level 0-1A
    if os.path.isfile(root + "/log/logfilePAR_01A.txt"):
        write_option = 'a'
    else:
        with open(root + "/log/logfilePAR_01A.txt", "w") as f:
            f.close()

    if infile in open(root + "/log/logfilePAR_01A.txt").read():
        next
    else:
        with open(root + "/log/logfilePAR_01A.txt", "a") as f:
            f.write(infile + "\n")
            f.close()

        # Initalize object with its attributes
        PAR = PAR_data()
        # Read data
        #infile = '../renku_ppmooring/ppmooring/data/Level0/RBR_PAR/RBR_PAR_1000cm/conv/PAR4-10m-079968_20201119.xlsx'

        PAR.read_data(infile)
        # Set output name and folder for Level 1A
        PAR.set_output(directories["Level1A_dir"], "L1A", "PAR")
        # First quality flag at level 0
        PAR.quality_flags_level0()
        # Store Level 1A data to NetCDF
        PAR.to_NetCDF("L1A")
        # Store Level 1A data to csv
        PAR.to_csv("L1A")

###Process level 1A - 1B
# concatenate files at level 1A and parse into 15 day period files

PAR_files.parsetimedata("Level1A", '15D')

# Initalize object with its attributes
PAR1B = PAR_data()

#Determine maintenance dates
##yesno = input("need to export maintenance.csv from PAR ? 1 yes, 0 no")
##if yesno == "1":

# Get list of files to read
dir_list = next(os.walk(directories["Level1A_dir"]))[1]
dir_list = ['LexplorePPMooringPAR_3000cm']

#Get maintenance dates from PAR 2000cm and PAR3000cm
df_final = pd.DataFrame()

for folder in dir_list:
    inpath = os.path.join(directories["Level1A_dir"],folder)
    filestoread = os.listdir(inpath)
    filestoread = fnmatch.filter(filestoread, '*.csv')

    vdepth = 3000

    valid_range = {
        "50": [0, 2000],
        "250": [0, 2000],
        "500": [0, 1500],
        "1000": [0, 1000],
        "2000": [0, 30],
        "3000": [0, 10]}

    dfall = pd.DataFrame()
    for file in filestoread:
        infile = os.path.join(inpath, file)
        dftemp = pd.read_csv(infile)

        # reset index
        dftemp = dftemp.reset_index(inplace=False, drop=True)
        time = dftemp.time
        PAR = dftemp.PAR

        # Check if variable is within valid range

        indlow = np.where(PAR >= valid_range[str(vdepth)][1])[0]
        indhigh = np.where(PAR >= valid_range[str(vdepth)][1])[0]

        timefirst = pd.DataFrame(time[indlow])
        timefirst = timefirst.rename(columns={'time':'debut'})
        timeend = pd.DataFrame(time[indhigh])
        timeend = timeend.rename(columns={'time': 'debut'})

        line = pd.concat([timefirst,timeend],axis=1)

        #dfall = dfall.append(line)
        dfall = pd.concat([dfall, line], ignore_index=True)

    df_final = pd.concat([df_final,dfall],ignore_index=True)

    # df_final = df_final.reset_index()
    # df_final["debut"] =pd.to_datetime(df_final["debut"],unit='s')
    # df_final["end"] = pd.to_datetime(df_final["end"], unit='s')

    filename = root + "/scripts/PAR_maintenance.csv"
    filename = filename.replace(os.sep, '/')
    df_final.to_csv(filename, sep=";", header=True, index=False)

#Reset filestoread
PAR1B.get_files_to_read_Leveltemp(directories["Level1A_dir"])

for infile in PAR1B.filestoread:
    #infile = '../renku_ppmooring/ppmooring/data/Level1A/LexplorePPMooringPAR_3000cm/temp/20190124_20190129_LexplorePPMooringPAR_3000cmLevel1A_temp.csv'
    infile = infile.replace(os.sep, '/')

    # Create log files containing treated files level 0-1A
    if os.path.isfile(root + "/log/logfilePAR_1A1B.txt"):
        write_option = 'a'
    else:
        with open(root + "/log/logfilePAR_1A1B.txt", "w") as f:
            f.close()

    if infile in open(root + "/log/logfilePAR_1A1B.txt").read():
        next
    else:
        # Read data and store it to NetCDF file
        # Initalize object with its attributes
        PAR = PAR_data()
        PAR.read_data_Leveln(infile)
        # Set output name and folder for Level 1B
        PAR.set_output(directories["Level1B_dir"], "L1B", "DO")

        ## PPMooring DO Data Level 0 --> Level 1A
        PAR.quality_flags_level0()

        ## PPMooring DO Data Level 1A --> Level 1B
        # Run quality assurance tests and flag each datapoint
        # Flagging bad data
        PAR.quality_flags_level1A()

        # Store Level 1B data to NetCDF
        PAR.to_NetCDF("L1A", "L1B")
        # Store Level 1B data to csv
        PAR.to_csv("L1A", "L1B")

        ## PPMooring DO Data Level 1b --> Level 2
        # Create Level 2 data, will be displayed on Datalakes

        # Mask bad data from Level 1B
        PAR.mask_data()
        # Store Level 2 data to NetCDF
        PAR.Level2_to_NetCDF(directories["Level2_dir"], "PAR")

        #Write on log files
        with open(root + "/log/logfilePAR_1A1B.txt", "a") as f:
            f.write(infile + "\n")
            f.close()


## DO Begin ######################################

# Initalize object with its attributes
DO_files = DO_data()
# Get list of files to read
DO_files.get_files_to_read_Level0(input_dirs["DO"], instruments = ["rbr", "exo", "minidot"])

##Check and delete empty files
for infile in DO_files.filestoread:
    fsize = os.stat(infile).st_size
    if fsize <= 1000: #delete files <1000 bytes
        os.remove(infile)

# Renew list of files to read
DO_files.get_files_to_read_Level0(input_dirs["DO"], instruments = ["rbr", "exo", "minidot"])


# Loop over all the files to read
for infile in DO_files.filestoread:
    ## PPMooring DO Data Level 0 --> Level 1A

    ##infile = DO_files.filestoread[1000]
    infile = infile.replace(os.sep,'/')

    # Create log files containing treated files level 0-1A
    if os.path.isfile(root + "/log/logfileDO_01A.txt"):
        write_option = 'a'
    else:
        with open(root + "/log/logfileDO_01A.txt", "w") as f:
            f.close()

    if infile in open(root + "/log/logfileDO_01A.txt").read():
        next
    else:
        with open(root + "/log/logfileDO_01A.txt", "a") as f:
            f.write(infile+ "\n")
            f.close()

        # Read data and store it to NetCDF file
        # Initalize object with its attributes
        DO = DO_data()

        # Read data
        if "rbr" in infile and "txt" in infile:
            DO.read_data_rbr(infile)
        elif "rbr" in infile and "xlsx" in infile:
            DO.read_data_rbr_excel(infile)
        elif "exo" in infile:
            DO.read_data_exo(infile)
        elif "minidot" in infile:
            DO.read_data_minidot(infile)

        # Set output name and folder for Level 1A
        DO.set_output(directories["Level1A_dir"], "L1A", "DO")

        # Store Level 1A data to NetCDF
        DO.to_NetCDF("L1A")

        # Store Level 1A data to csv
        DO.to_csv("L1A")


###Process level 1A - 1B
# concatenate files at level 1A and parse into 15 day period files
yesno = input("need to export to 15days dataframes for DO data? 1 yes, 0 no")
if yesno =="1":
    DO_files.parsetimedata("Level1A",'15D')

# Initalize object with its attributes
DO1B = DO_data()

# Get list of files to read
DO1B.get_files_to_read_Leveltemp(directories["Level1A_dir"])

for infile in DO1B.filestoread:
    infile = infile.replace(os.sep, '/')

    # Create log files containing treated files level 0-1A
    if os.path.isfile(root + "/log/logfileDO_1A1B.txt"):
        write_option = 'a'
    else:
        with open(root + "/log/logfileDO_1A1B.txt", "w") as f:
            f.close()

    if infile in open(root + "/log/logfileDO_1A1B.txt").read():
        next
    else:
        # Read data and store it to NetCDF file
        # Initalize object with its attributes
        DO = DO_data()
        DO.read_data_Leveln(infile)
        # Set output name and folder for Level 1B
        DO.set_output(directories["Level1B_dir"], "L1B", "DO")

        ## PPMooring DO Data Level 0 --> Level 1A
        DO.maintenance = DO_data.DO_maintenance_extract()
        DO.quality_flags_level0()

        ## PPMooring DO Data Level 1A --> Level 1B
        # Run quality assurance tests and flag each datapoint
        # Flagging bad data
        DO.quality_flags_level1A()

        # Store Level 1B data to NetCDF
        DO.to_NetCDF("L1A", "L1B")
        # Store Level 1B data to csv
        DO.to_csv("L1A", "L1B")

        ## PPMooring DO Data Level 1b --> Level 2
        # Create Level 2 data, will be displayed on Datalakes

        # Mask bad data from Level 1B
        DO.mask_data()
        # Store Level 2 data to NetCDF
        DO.Level2_to_NetCDF(directories["Level2_dir"], "DO")

        #Write on log files
        with open(root + "/log/logfileDO_1A1B.txt", "a") as f:
           f.write(infile + "\n")
           f.close()
## DO End ######################################


## Temperature begin ##########################

# Initalize object with its attributes
Temp_files = Temperature_data()
# Get list of files to read
Temp_files.get_files_to_read_Level0(input_dirs["Temp"])
# Loop over all the .txt files to read
for infile in Temp_files.filestoread:

    infile = infile.replace(os.sep,'/')

    # Create log files containing treated files level 0-1A
    if os.path.isfile(root + "/log/logfileTemp_01A.txt"):
        write_option = 'a'
    else:
        with open(root + "/log/logfileTemp_01A.txt", "w") as f:
            f.close()

    if infile in open(root + "/log/logfileTemp_01A.txt").read():
        next
    else:
        with open(root + "/log/logfileTemp_01A.txt", "a") as f:
            f.write(infile + "\n")
            f.close()

        # Read data and store it to NetCDF file

        # Initalize object with its attributes
        Temp = Temperature_data()
        # Read data
        if ".txt" in infile:
            Temp.read_data_txt(infile)
        elif ".csv" in infile:
            Temp.read_data_csv(infile)
        # Set output name and folder for Level 1A
        Temp.set_output(directories["Level1A_dir"], "L1A", "Temp")
        # Store Level 1A data to NetCDF
        Temp.to_NetCDF("L1A")

        # Store Level 1A data to csv
        Temp.to_csv("L1A")

###Process level 1A - 1B
# concatenate files at level 1A and parse into 15 day period files
yesno = input("need to export to 15days dataframes for Temp data? 1 yes, 0 no")
if yesno =="1":
    Temp_files.parsetimedata("Level1A",'15D')

# Initalize object with its attributes
Temp1B = Temperature_data()

# Get list of files to read
Temp1B.get_files_to_read_Leveltemp(directories["Level1A_dir"])

for infile in Temp1B.filestoread:
    infile = infile.replace(os.sep, '/')

    # Create log files containing treated files level 0-1A
    if os.path.isfile(root + "/log/logfileTemp_1A1B.txt"):
        write_option = 'a'
    else:
        with open(root + "/log/logfileTemp_1A1B.txt", "w") as f:
            f.close()

    if infile in open(root + "/log/logfileTemp_1A1B.txt").read():
        next
    else:
        # Read data and store it to NetCDF file
        # Initalize object with its attributes
        Temp = Temperature_data()
        Temp.read_data_Leveln(infile)
        # Set output name and folder for Level 1B
        Temp.set_output(directories["Level1B_dir"], "L1B", "DO")

        ## PPMooring DO Data Level 0 --> Level 1A
        Temp.maintenance = Temperature_data.Temp_maintenance_extract()
        Temp.quality_flags_level0()

        ## PPMooring DO Data Level 1A --> Level 1B
        # Run quality assurance tests and flag each datapoint
        # Flagging bad data
        Temp.quality_flags_level1A()

        # Store Level 1B data to NetCDF
        Temp.to_NetCDF("L1A", "L1B")
        # Store Level 1B data to csv
        Temp.to_csv("L1A", "L1B")

        ## PPMooring DO Data Level 1b --> Level 2
        # Create Level 2 data, will be displayed on Datalakes

        # Mask bad data from Level 1B
        Temp.mask_data()
        # Store Level 2 data to NetCDF
        Temp.Level2_to_NetCDF(directories["Level2_dir"], "DO")

        #Write on log files
        with open(root + "/log/logfileTemp_1A1B.txt", "a") as f:
           f.write(infile + "\n")
           f.close()

## Temperature end ##########################
print("ok")