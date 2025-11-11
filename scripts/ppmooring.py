# -*- coding: utf-8 -*-

#### Import Packages
import re

import netCDF4
import pandas as pd
import numpy as np
import fnmatch
import os
import datetime as dt
from datetime import datetime, timezone
from copy import deepcopy
import matplotlib.pyplot as plt
from collections import defaultdict
from functools import partial
from tqdm import tqdm
from tsmoothie.smoother import *
from sklearn.cluster import KMeans


#### Define Classes
class ppmooring(object):
    ppmooringroot = '../renku_ppmooring/ppmooring'
    ## Parent class for all the different sensors of the ppmooring module

    def set_output(self, output_dir, level, instrument):
        ## Set output name and create output folder it not already existing
       
        # Set output folder and output file name
        if level == "L2":    
            outfolder = output_dir + self.Level2_outfolder
            date_time = datetime.utcfromtimestamp(self.time_grid.iloc[0]).strftime('%Y%m%d_%H%M%S')
            self.outfile = outfolder + "/" + level + "_" + instrument + "_" + date_time + ".nc"
        elif level == "L1B":
            outfolder = output_dir + self.outfolders[self.folder]

            preoutfile = self.name1B.replace("temp/","")
            preoutfile = preoutfile.replace("_temp", "")
            preoutfile = preoutfile.replace("1A", "1B")
            preoutfile = preoutfile[0:-4]

            self.outfile = preoutfile + '.nc'
            self.outfilecsv = preoutfile + '.csv'
        else:
            outfolder = output_dir + self.outfolders[self.folder]
            date_time = datetime.utcfromtimestamp(self.time.iloc[0]).strftime('%Y%m%d_%H%M%S')
            self.outfile = outfolder + "/" + level + "_" + instrument + "_" + \
                str(self.depth) + "cm" + "_" + date_time + ".nc"
            self.outfilecsv = outfolder + "/" + level + "_" + instrument + "_" + \
                              str(self.depth) + "cm" + "_" + date_time + ".csv"
        
        # Make sure the output folder exists
        if not os.path.exists(outfolder):
            os.makedirs(outfolder)
            
    def to_NetCDF(self, *args):
        ## Save data to NetCDF file
        
        # Check which variables to store in NetCDF file
        if args:
            self.var_dict = {}
        else:
            print("No variables will be stored to NetCDF file")
            
        # Store variables and attributes
        for level in args:
            self.var_dict.update(getattr(self, level+"_dict"))
                
        # Sort dict alphabetically
        sorted_dict = {}
        for key, values in sorted(self.var_dict.items()): 
            sorted_dict[key] = values 
        self.var_dict = sorted_dict
            
        # Make sure NetCDF files are closed
        try: ncfile.close() 
        except: pass
                    
        # Check if file already exists
        if os.path.isfile(self.outfile):
            # Open file in 'append' mode
            ncfile = netCDF4.Dataset(self.outfile,mode='a',format='NETCDF4_CLASSIC')
            
            # Add Level 2 data
            if "L2" in args:
                # Add data to variables
                for key, values in self.var_dict.items():
                    self.addDataNetCDF(getattr(self,key),0,ncfile,"L2",**values)
                    
                # Close NetCDF file
                ncfile.close()
            
            # Add Level 1A/1B data
            else:
                # Check if last timestamp has already been added, if not data will be added
                if self.time.iloc[-1] in ncfile.variables['time'][:]:
                    print("Duplicated run, no data added")
                    ncfile.close() 
                    return
                
                else:
                    # Get index of last saved timestamp
                    time = ncfile.variables['time']
                    idx = len(time[:])
                    
                    # Add data to coordinate variable 'time'
                    time[:] = np.append(time[:], self.time)
                    
                    # Add data to variables
                    for key, values in self.var_dict.items():
                        self.addDataNetCDF(getattr(self,key),idx,ncfile,**values)
                        
                    # Close NetCDF file
                    ncfile.close()
                
        else:           
            # Create new NetCDF file
            ncfile = netCDF4.Dataset(self.outfile, mode='w',format='NETCDF4_CLASSIC') 
            
            # Add global attributes
            for key in self.att:
                setattr(ncfile, key, self.att[key])
    
            # Create Dimensions
            for key, values in self.dim_dict.items():
                ncfile.createDimension(values['dim_name'], values['dim_size'])
                
            # Add data to coordinate variables
            for key, values in self.dim_var_dict.items():
                self.addVarNetCDF(getattr(self,key),ncfile,**values)
                
            # Add data to variables
            if "L2" in args:
                # Create variable for each depth
                for var in [*self.L2_dict]:
                    temp_dict = deepcopy(self.L2_dict[var])
                    for entry in list(self.depthdict.values()):
                        temp_dict["var_name"] = str(entry/100) + "m"
                        
                        var_temp = ncfile.createVariable(temp_dict["var_name"], np.float64, temp_dict["dim"])
                        var_temp.units = temp_dict["unit"]
                        var_temp.long_name = temp_dict["longname"]
                
                # Add data of current depth
                for key, values in self.var_dict.items():
                    self.addDataNetCDF(getattr(self,key),0,ncfile,"L2",**values)    
            else:
                for key, values in self.var_dict.items():
                    self.addVarNetCDF(getattr(self,key),ncfile,**values)
                
            # Close NetCDF file
            ncfile.close()
    
    def addDataNetCDF(self,data,idx,ncfile,*args,**kwargs):
        ## Add data to existing NetCDF file (only along time dimension)
         
        # Check if variable exists, if not create new one
        try: 
            var = ncfile.variables[kwargs["var_name"]]
        except:
            var = ncfile.createVariable(kwargs["var_name"], np.float64, kwargs["dim"])
            var.units = kwargs["unit"]
            var.long_name = kwargs["longname"]
        
        # Add Level 2 data if args specified, otherwise Level 1A/1B
        if args:   
            # Add data to variable
            var[self.idx_grid] = data
            
        else:
            # Check dimension size
            if "dim" in kwargs.keys():
                ndim = len(kwargs["dim"])
            else:
                ndim = 0        
                
            # Add data to variable
            if ndim == 2:
                var[:,idx:] = data
            elif ndim == 1:
                var[idx:] = data
            else:
                var[:] = data
            
    def addVarNetCDF(self,data,ncfile,**kwargs):
        ## Create variable in new NetCDF file
        
        # Create variable and check for dimension variable
        if "dim" in kwargs.keys():
            var = ncfile.createVariable(kwargs["var_name"], np.float64, kwargs["dim"])
            opt = 1
        else:
            var = ncfile.createVariable(kwargs["var_name"], np.float64)
            opt = 0
           
        # Add unit and longname
        var.units = kwargs["unit"]
        var.long_name = kwargs["longname"]
               
        # Add data to variable
        if opt ==1:
            var[:] = data # vector with dimension >1
        elif isinstance(data,int) == True:
            var[:] = data # depth at level 0
        else:
            var[:] = data[0] #depth at level 1A temp
        
    def mask_data(self):
        ## Replace bad values from Level 1B by NaN
        
        for var in [*self.L2_dict]:
            idx = getattr(self, var+"_qual") > 0
            
            getattr(self, var)[idx] = np.nan
            
    def grid_data(self, start_period, end_period):
        ## Store data in grid before outputting it as Level 2
        
        # Define time grid
        time_grid = np.array([start_period])
        while (time_grid[-1] + dt.timedelta(minutes = 5)) < end_period:
            time_grid = np.append(time_grid, time_grid[-1] + dt.timedelta(minutes = 5))  
         
        # Get the days corresponding to the timestamps of the measurements
        date_time = np.array([datetime.utcfromtimestamp(time) for time in self.time])
        
        # Consider duplicate time values because of change from winter time to summer time
        date_time, idx_unique = np.unique(date_time, return_index=True)
        for var in [*self.L2_dict]:
            setattr(self, var, getattr(self, var).iloc[idx_unique])
            setattr(self, var, getattr(self, var).reset_index(drop = True))
            
        # Check if measurments exist within 5 minutes window of any time_grid value
        self.idx_grid = np.full(time_grid.shape, False)
        idx_data = np.full(date_time.shape, False)
        for idx, entry in enumerate(date_time):
               
            # Find closest time_grid value for current entry
            timedelta = abs(time_grid - entry)
            idx_min = np.argmin(timedelta)
            min_timedelta = timedelta[idx_min]
            
            # Check if there's a matching timestamp
            if min_timedelta == dt.timedelta():
                
                self.idx_grid[idx_min] = True
                idx_data[idx] = True
                
            # Check if timestamp within 5 minutes window corresponds to minimum
            elif min_timedelta < dt.timedelta(seconds = 150):
                  
                candidat = time_grid[idx_min]
                if idx == 0:
                    if min_timedelta < abs(candidat - date_time[idx+1]):
                    
                        self.idx_grid[idx_min] = True
                        idx_data[idx] = True
                        
                elif 0 < idx and idx < len(date_time) - 1:
                    if min_timedelta < abs(candidat - date_time[idx-1]) and min_timedelta < abs(candidat - date_time[idx+1]):
                
                        self.idx_grid[idx_min] = True
                        idx_data[idx] = True
                
                elif idx == len(date_time) - 1:
                    if min_timedelta < abs(candidat - date_time[idx-1]):
                    
                        self.idx_grid[idx_min] = True
                        idx_data[idx] = True
                
        # Update time and measurement data
        self.time_grid = pd.Series([value.replace(tzinfo=timezone.utc).timestamp() for value in time_grid])
        for var in [*self.L2_dict]:
            setattr(self, var, getattr(self, var)[idx_data])
            setattr(self, var, getattr(self, var).reset_index(drop = True))
        
    def Level2_to_NetCDF(self, directory, instrument):
        ## Save Level2 data to NetCDF files where the data is seperated
        ## s.t. each file spans over a time period of 10 days
        
        # Define time interval in days 
        interval = 10
        
        # Get the datetime objects corresponding to the timestamps of the measurements
        date_time = np.array([datetime.utcfromtimestamp(time) for time in self.time])
        
        # Set beginning of time grid
        grid_begin = datetime(2018, 1, 1)
                    
        # Set end of time grid
        grid_end = date_time[-1]
        
        # Define time grid
        time_grid = np.array([grid_begin])
        while time_grid[-1] <= grid_end:
            time_grid = np.append(time_grid, time_grid[-1] + dt.timedelta(days = interval))
        
        # Get index of measurements for each period
        idx_period = np.arange(1, len(time_grid))     
        for period in idx_period:
            
            # Slice values to store
            idx = np.where((time_grid[period-1] <= date_time) & (date_time < time_grid[period]))[0]
            
            # Check if any measurements in current period
            if len(idx) == 0:
                continue
            else:   
                # Create copy of current class instance to save current period to NetCDF
                period_obj = deepcopy(self)
                
                # Slice data
                period_obj.time = period_obj.time[idx]
                for var in [*period_obj.L2_dict]:
                    setattr(period_obj, var, getattr(period_obj,var)[idx])
                    period_obj.L2_dict[var]["var_name"] = str(period_obj.depth[0]/100) + "m"
                    
                # Grid data
                period_obj.grid_data(time_grid[period-1], time_grid[period])
                    
                # Change dictonary key to guarantee compability with the .to_NetCDF method
                period_obj.dim_dict["time_grid"] = period_obj.dim_dict.pop("time")
                period_obj.dim_var_dict["time_grid"] = period_obj.dim_var_dict.pop("time")
                
                # Set output name
                period_obj.set_output(directory, "L2", instrument)            
                
                # Save to NetCDF
                period_obj.to_NetCDF("L2")

    def to_csv(self,*args):
        ## Save data to csv files
        # Check which variables to store in csv file
        if args:
            self.var_dict = {}
        else:
            print("No variables will be stored to CSV file")

        # Store variables and attributes
        for level in args:
            self.var_dict.update(getattr(self, level + "_dict"))

        # Sort dict alphabetically
        sorted_dict = {}
        for key, values in sorted(self.var_dict.items()):
            sorted_dict[key] = values
        self.var_dict = sorted_dict

        # Check if file already exists
        if os.path.isfile(self.outfilecsv):
            print("Duplicated run, no csv created")
        else:
            data = {'time' : getattr(self,'time')}
            for key,values in self.var_dict.items():
                data[key] = getattr(self,key)

            data = pd.DataFrame.from_dict(data)
            data.to_csv(self.outfilecsv,index=False)

class DO_data(ppmooring):
    ## Class for the DO data

    def __init__(self):
        ## Set attributes for this class

        #root
        self.root = ppmooring.ppmooringroot
        # Output folders
        self.outfolders = {
            "DO_50cm": "/LexplorePPMooringDO_50cm",
            "DO_250cm": "/LexplorePPMooringDO_250cm",
            "DO_500cm": "/LexplorePPMooringDO_500cm",
            "DO_1000cm": "/LexplorePPMooringDO_1000cm",
            "DO_1500cm": "/LexplorePPMooringDO_1500cm",
            "DO_2000cm": "/LexplorePPMooringDO_2000cm",
            "DO_3000cm": "/LexplorePPMooringDO_3000cm",
            "DO_5000cm": "/LexplorePPMooringDO_5000cm",
            "DO_10000cm": "/LexplorePPMooringDO_10000cm"
        }

        self.Level2_outfolder = "/LexplorePPMooringDO"
        
        # Matrix depth - folder
        self.depthdict = {
            "DO_50cm": 50,
            "DO_250cm": 250,
            "DO_500cm": 500,
            "DO_1000cm": 1000,
            "DO_1500cm": 1500,
            "DO_2000cm": 2000,
            "DO_3000cm": 3000,
            "DO_5000cm": 5000,
            "DO_10000cm": 10000
        }
        
        # General attributes
        self.att = {
            "institution": "EPFL",
            "source": "Mooring DO",
            "references": "LéXPLORE commun instruments Viet Tran Khac <viet.tran-khac@inra.fr>",
            "history": "See history on Renku",
            "conventions": "CF 1.7",
            "comment": "Dissolved Oxygen from Mooring M2 on Lexplore Platform in Lake Geneva",
            "title": "Lexplore Mooring DO"
        }
        
        # Dimension attributes
        self.dim_dict = {
            'time': {'dim_name':'time', 'dim_size':None}
        }
        
        # Coordinate variable attributes
        self.dim_var_dict = {
            'time': {'var_name':'time', 'dim':('time',), 'unit':'seconds since 1970-01-01 00:00:00', 'longname':'time'}
        }
        
        # Level 1A - Variable attributes
        self.L1A_dict = {
            'depth': {'var_name':'depth','unit':'cm', 'longname':'depth'},
            'O2Sat': {'var_name':'O2Sat', 'dim':('time',), 'unit':'%', 'longname':'oxygen saturation'},
            'O2Sat_qual': {'var_name': 'O2Sat_qual', 'dim': ('time',),
                           'unit': '0 = nothing to report, 1 = more investigation',
                           'longname': 'oxygen saturation quality flag'},
            'O2': {'var_name':'O2', 'dim':('time',), 'unit':'mgO2/L', 'longname':'dissolved_oxygen'},
            'O2_qual': {'var_name': 'O2_qual', 'dim': ('time',),
                        'unit': '0 = nothing to report, 1 = more investigation',
                        'longname': 'dissolved_oxygen quality flag'}
        }
        
        # Level 1B - Qualtiy flag attributes
        self.L1B_dict = {
            'O2Sat_qual': {'var_name':'O2Sat_qual', 'dim':('time',), 'unit':'0 = nothing to report, 1 = more investigation', 'longname':'oxygen saturation quality flag'},
            'O2_qual': {'var_name':'O2_qual', 'dim':('time',), 'unit':'0 = nothing to report, 1 = more investigation', 'longname':'dissolved_oxygen quality flag'}
        }
        
        # Level 2 - Variable attributes
        self.L2_dict = {
            'O2': {'var_name':'O2', 'dim':('time',), 'unit':'mgO2/L', 'longname':'dissolved_oxygen'}
        }
     
    def get_files_to_read_Level0(self, input_dir, instruments):
        ## Get a list of all files to read

        # List all subfolders
        folders = os.listdir(input_dir)
        folders = fnmatch.filter(folders, 'DO_*')
        folders = [input_dir + "/" + s for s in folders]
        
        self.filestoread = []
        # Loop through all subfolders
        for idxfolder, folder in enumerate(folders):
            subfolders = os.listdir(folder)
            
            for instrument in instruments:  
                if instrument in subfolders:
                    files_subfolder = [folder + "/" + instrument + "/" + s for s in os.listdir(folder + "/" + instrument)]
                    self.filestoread = np.append(self.filestoread, files_subfolder)

    def get_files_to_read_Leveltemp(self, input_dir):
        ## Get a list of all files to read at level n

        # List all subfolders
        folders = os.listdir(input_dir)
        folders = fnmatch.filter(folders, 'LexplorePPMooringDO*')
        folders1 = [input_dir + "/" + s + "/temp" for s in folders]

        self.filestoread = []
        # Loop through all subfolders
        for folder in (folders1):
            files_subfolder = os.listdir(folder)
            files_subfolder = [folder + "/" + s for s in os.listdir(folder)]
            self.filestoread = np.append(self.filestoread, files_subfolder)

    def read_data_rbr(self, infile):
        ## Read RBR data
        
        # Store infile variable
        self.infile = infile
        
        # Get input folder
        self.folder = self.infile.rsplit("/")[-3]
                
        # Import metadata and get indicies of metadata
        with open(self.infile, encoding="utf8", errors='ignore') as f:
            for idx, line in enumerate(f):
                if 'Serial' in line:
                    idx_serial = idx
                if 'Channel[1].calibration' in line:
                    idx_calibr = idx
                if 'NumberOfSamples' in line:
                    idx_meta = idx + 1
                    break
            f.seek(0)
            lines = f.readlines()
            line_name = lines[idx_meta + 1]
            meta_data = lines[0:idx_meta]
        
        # Extract meta-data
        numsensor = meta_data[idx_serial].rsplit("=")[1].rstrip("\n")
        numcalibration = meta_data[idx_calibr].rsplit("=")[1].rstrip("\n")
    
        # Read data
        #20221201-insert encoding = 'unicode_escape', because utf-8 encoding bugs

        df = pd.read_csv(self.infile, delim_whitespace=True, header=None, skiprows=idx_meta + 2, encoding='unicode_escape')
        
        # Take into account that files do not always have the same columns
        if "doxy07" in line_name or "doxy01" in line_name:
            df = df.drop([3, 5, 6], axis=1)
            df.columns = ["Date", "Time", "Temp", "O2Sat", "O2"]
            
            cols = ['O2Sat', 'O2']
            df[cols] = df[cols].apply(pd.to_numeric, errors='coerce', axis=1)
        else:
            df = df.drop([4], axis = 1)
            df.columns = ["Date", "Time", "O2", "Temp", "O2Sat"]
            
            cols = ['O2Sat', 'O2']
            df[cols] = df[cols].apply(pd.to_numeric, errors='coerce', axis=1)
            
            # Convert O1 µmol/L to mg/L&4
            df["O2"] = df["O2"] * 32 / 1000
            
        df = df.reset_index(drop=True)
    
        # Replace date and time by timestamp
        datehour = df["Date"] + " " + df["Time"]
        datehour = datehour.apply(lambda x: datetime.strptime(x, '%Y-%m-%d %H:%M:%S.%f'))
        df["Timestamp"] = [datetime.timestamp(s) for s in datehour]
    
        df = df[["Timestamp", "O2Sat", "O2"]]
        
        # Sort by date and reset index
        df = df.sort_values(by=['Timestamp'])
        self.df = df.reset_index(inplace=False, drop=True)
        
        # Sensor data, depth and point
        self.meta_data = {
            "sensor": numsensor,
            "calibration": numcalibration,
            "instrument": "rbr",
            "depth": self.depthdict[self.folder],
            "point": self.folder
        }

        # Store variables
        self.time = df["Timestamp"] # need to be checked      
        self.depth = self.meta_data["depth"]
        self.O2Sat = df["O2Sat"]
        self.O2 = df["O2"]

        # Create vector quality for data quality flags at Level 1A
        self.O2_qual = np.zeros(len(self.O2))
        self.O2Sat_qual = np.zeros(len(self.O2_qual))

    def read_data_rbr_excel(self, infile):
        ## Read RBR data
        # Store infile variable
        self.infile = infile

        # Get input folder
        self.folder = self.infile.rsplit("/")[-3]

        # Import metadata and get indicies of metadata
        metadata = pd.read_excel(open(infile, 'rb'), sheet_name='Metadata')

        idx = np.where(metadata.iloc[:, 0].str.contains("Serial") == True)[0][0] + 1
        numsensor = metadata.iloc[idx, 0]

        idx = np.where(metadata.iloc[:, 0].str.contains("Atmospheric pressure") == True)[0][0] + 1
        atmospheric_pressure = metadata.iloc[idx, 0]
        density = metadata.iloc[idx, 1]
        meta_Temperature = metadata.iloc[idx, 3]
        pressure = metadata.iloc[idx, 4]
        meta_conductivity = metadata.iloc[idx, 5]
        salinity = metadata.iloc[idx, 6]
        conductivity_specific = metadata.iloc[idx, 7]

        # Read data
        df = pd.read_excel(open(infile, 'rb'), sheet_name='Data', skiprows=1)

        ncol = df.shape[1]
        # Take into account that files do not always have the same columns
        if ncol ==7:
            #Data rbr with Temperature	Pressure	O₂ saturation	Sea pressure	Depth	 O₂ concentration (umol/L)
            df.columns = ["Date","Temp","Pressure","O2Sat","Sea_pressure","Depth","O2"]

            ##df["O2"] =(df["O2"] +30.06) * 0.032 #Convert to mg/L
            termA = np.exp(7.7117 - 1.31403 * np.log((df["Temp"] + 45.93)))
            termB = np.exp(11.8571-(3840.7 / (df["Temp"]+273.15))-(216961/np.power((df["Temp"] + 273.15),2)))
            termC = (0.000975 - (0.00001426 * df["Temp"]) + (0.00000006436 * np.power(df["Temp"],2)))
            termABC = (termA * (df["Pressure"] / 10) * (1 - termB / df["Pressure"] / 10) * (1 - termC * df["Pressure"]))
            cpO2 = (termABC / (1 - np.exp(11.8571 - (3840.7 / (df["Temp"] + 273.15)) - (216961 / np.power((df["Temp"] + 273.15),2)))) / (1 - (0.000975 - (0.00001426 * df["Temp"]) + (0.00000006436 * np.power(df["Temp"],2)))))

            df["O2"]= df["O2Sat"]/100 * cpO2

            #remove all rows < or > than depth
            vdepth = self.depthdict[self.folder]
            idx = np.where((df["Depth"]<vdepth/100*0.8) | (df["Depth"]>vdepth/100*1.2))[0]
            df["O2"][idx]=0
            df["O2Sat"][idx]=0


            df = df.reset_index(drop=True)
            df["Timestamp"] = [datetime.timestamp(s) for s in df['Date']]

        if ncol == 5:
            df.columns = ["Date","O2","Temp","Corrected_phase","O2Sat"]
            df["O2"] = (df["O2"]+df["Corrected_phase"])*0.032
            ##df["O2Sat"] = df["O2Sat"]*1.0502 #1/p correction 400m

            df = df.reset_index(drop=True)

            # Replace date and time by timestamp
            df["Timestamp"] = [datetime.timestamp(s) for s in df['Date']]

        df = df[["Timestamp", "O2Sat", "O2"]]

        # Sort by date and reset index
        df = df.sort_values(by=['Timestamp'])
        self.df = df.reset_index(inplace=False, drop=True)

        # Sensor data, depth and point
        self.meta_data = {
            "sensor": numsensor,
            "atmospheric_pressure": atmospheric_pressure,
            "density": density,
            "instrument": "rbr",
            "meta_conductivity": meta_conductivity,
            "depth": self.depthdict[self.folder],
            "point": self.folder
        }

        # Store variables
        self.time = df["Timestamp"]  # need to be checked
        self.depth = self.meta_data["depth"]
        self.O2Sat = df["O2Sat"]
        self.O2 = df["O2"]

        # Create vector quality for data quality flags at Level 1A
        self.O2_qual = np.zeros(len(self.O2))
        self.O2Sat_qual = np.zeros(len(self.O2_qual))

    def read_data_exo(self, infile):
        ## Read EXO data    
            
        # Store infile variable
        self.infile = infile
        
        # Get input folder
        self.folder = self.infile.rsplit("/")[-3]
            
        # Import metadata and get indicies of metadata
        with open(self.infile, encoding="utf8", errors='ignore') as f:
            for idx, line in enumerate(f):
                if 'Sonde ID' in line:
                    idx_sonde = idx
                    break
            f.seek(0)
            lines = f.readlines()
            meta_data = lines[0:26]

        # Extract meta-data
        numsensor = meta_data[idx_sonde].rsplit(";")[1].rsplit(" ")[1]    
        numcalibration = -99
        
        # Read data
        df = pd.read_csv(self.infile, sep=";", skiprows=26, header=None, encoding='unicode_escape')
        
        ncol = len(df.columns)
        # If data has 26 columns = no TSS, if 27 columns = TSS column included
        if ncol == 27:
            df = df.drop([2, 3, 4, 6], axis=1)
            df.columns = ["Date", "Time", "Batt", "pH", "pHmV", "NNO3mv", "NNO3",
                          "ChlaRFU", "Chla", "PCRFU", "PC", "Temp", "Cond", "SpCond", "Sal", "nLFCond",
                          "TDS", "Turbidity", "TSS", "O2Sat", "O2", "Pressure", "Depth"]
            
            # Add timestamp
            datehour = df["Date"] + " " + df["Time"]
            datehour = datehour.apply(lambda x: datetime.strptime(x, '%m/%d/%Y %H:%M:%S'))
            df["Timestamp"] = [datetime.timestamp(s) for s in datehour]
    
            df = df[["Timestamp", "Batt", "pH", "pHmV", "NNO3mv", "NNO3",
                     "ChlaRFU", "Chla", "PCRFU", "PC", "Temp", "Cond", "SpCond", "Sal", "nLFCond",
                     "TDS", "Turbidity", "TSS", "O2Sat", "O2", "Pressure", "Depth"]]
        else:
            df = df.drop([2, 3, 4, 6], axis=1)
            df.columns = ["Date", "Time", "Batt", "pH", "pHmV", "NNO3mv", "NNO3",
                          "ChlaRFU", "Chla", "PCRFU", "PC", "Temp", "Cond", "SpCond", "Sal", "nLFCond",
                          "TDS", "Turbidity", "O2Sat", "O2", "Pressure", "Depth"]
            df["TSS"] = -99
            
            # Add timestamp
            datehour = df["Date"] + " " + df["Time"]
            datehour = datehour.apply(lambda x: datetime.strptime(x, '%m/%d/%Y %H:%M:%S'))
            df["Timestamp"] = [datetime.timestamp(s) for s in datehour]
    
            df = df[["Timestamp", "Batt", "pH", "pHmV", "NNO3mv", "NNO3",
                     "ChlaRFU", "Chla", "PCRFU", "PC", "Temp", "Cond", "SpCond", "Sal", "nLFCond",
                     "TDS", "Turbidity", "TSS", "O2Sat", "O2", "Pressure", "Depth"]]
        
        df =df[["Timestamp", "O2Sat", "O2"]]
        
        # Sort by date and reset index
        df = df.sort_values(by=['Timestamp'])
        self.df = df.reset_index(inplace=False, drop=True)
        
        # Sensor data, depth and point
        self.meta_data = {
            "sensor": numsensor,
            "calibration": numcalibration,
            "instrument": "exo",
            "depth": self.depthdict[self.folder],
            "point": self.folder
        }

        # Store variables
        self.time = df["Timestamp"] # need to be checked      
        self.depth = self.meta_data["depth"]
        self.O2Sat = df["O2Sat"]
        self.O2 = df["O2"]

        # Create vector quality for data quality flags at Level 1A
        self.O2_qual = np.zeros(len(self.O2))
        self.O2Sat_qual = np.zeros(len(self.O2_qual))
        
    def read_data_minidot(self, infile):
        ## Read Minidot data       
        
        # Store infile variable
        self.infile = infile
        
        # Get input folder
        self.folder = self.infile.rsplit("/")[-3]
               
        # Import metadata and get indicies of metadata
        with open(self.infile, encoding="utf8", errors='ignore') as f:
            lines = f.readlines()
            meta_data = lines[0:3]
            
        # Extract meta-data
        numsensor = meta_data[0].rsplit("\n")[0]
        numcalibration = meta_data[1].rsplit(":")[2].rsplit('\n')[0]

        # Read data
        df = pd.read_csv(self.infile, sep=",", skiprows=3, header=None)
        df.columns= ["Timestamp", "Batt", "Temp", "O2", "Ratio"]
        
        # O2Sat column
        #Calculer O2Sat data
        DO_cp = 14.64 - (0.4227 * df["O2"]) + (0.009937 * np.power(df["O2"],2)) - (0.0001575 * np.power(df["O2"],3)) + (0.000001125 * np.power(df["O2"], 4))
        df["O2Sat"] = (df["O2"]/DO_cp)*100

        df =df[["Timestamp", "O2Sat", "O2"]]
        
        # Sort by date and reset index
        df = df.sort_values(by=['Timestamp'])
        self.df = df.reset_index(inplace=False, drop=True)

        # Sensor data, depth and point
        self.meta_data = {
            "sensor": numsensor,
            "calibration": numcalibration,
            "instrument": "minidot",
            "depth": self.depthdict[self.folder],
            "point": self.folder
        }

        # Store variables
        self.time = df["Timestamp"] # need to be checked      
        self.depth = self.meta_data["depth"]
        self.O2Sat = df["O2Sat"]
        self.O2 = df["O2"]

        # Create vector quality for data quality flags at Level 1A
        self.O2_qual = np.zeros(len(self.O2))
        self.O2Sat_qual = np.zeros(len(self.O2_qual))

    def read_data_Leveln(self,infile):
        # Read data
        df = pd.read_csv(infile,sep=";")# delim_whitespace=True, header=None, skiprows=idx_meta + 2)
        df = df.reset_index(drop=True)

        # Replace date and time by timestamp
        df["Timestamp"] = df["time"]

        df = df[["Timestamp", "O2Sat", "O2","O2Sat_qual", "O2_qual","depth"]]

        # Sort by date and reset index
        df = df.sort_values(by=['Timestamp'])
        self.df = df.reset_index(inplace=False, drop=True)

        # Store variables
        # Get input folder
        self.name1B = infile
        self.folder = "DO_"+str(df["depth"][0])+"cm"

        self.time = df["Timestamp"]  # need to be checked
        self.depth = df["depth"]
        self.O2Sat = df["O2Sat"]
        self.O2 = df["O2"]
        self.O2_qual = df["O2_qual"]
        self.O2Sat_qual = df["O2Sat_qual"]

    def parsetimedata(self, atlevel, freqday):
        root = self.root
        input_dir = os.path.join(root, "data", atlevel)
        input_dirs = os.listdir(input_dir)
        input_dirs = fnmatch.filter(input_dirs, 'LexplorePPMooringDO_*')

        for folder in input_dirs:
            depthfolder = os.path.join(input_dir, folder)
            files = os.listdir(depthfolder)
            files = files[0:(len(files) - 1)]  # remove temp folder
            files = fnmatch.filter(files,'*.csv')

            dfall = pd.DataFrame()
            for file in files:
                filepath = os.path.join(depthfolder, file)
                df = pd.read_csv(filepath)
                # dfall = dfall_append(df) #pandas <2.0
                dfall = pd.concat([dfall, df], ignore_index=True)

            # Split into 15 days files
            dfall['time1'] = dfall.time.apply(lambda x: datetime.fromtimestamp(x))

            # Split df_final into 15 day dataframes and expor temporarily into csv
            # freqday = '15D'
            days = [g for n, g in dfall.groupby(pd.Grouper(key="time1", freq=freqday))]

            for day in range(0, len(days)):
                if len(days[day]) == 0:
                    next
                else:
                    data = days[day]
                    data = data.sort_values(by=['time'])
                    data = data.reset_index(inplace=False, drop=True)

                    outfolder = os.path.join(depthfolder, "temp")
                    if not os.path.exists(outfolder):
                        os.makedirs(outfolder)

                    prefix_start = datetime.fromtimestamp(data["time"].iloc[0]).strftime('%Y%m%d')
                    prefix_end = datetime.fromtimestamp(data["time"].iloc[-1]).strftime('%Y%m%d')

                    datename = prefix_start + "_" + prefix_end

                    filename = outfolder + "/" + datename + "_" + folder + "Level1A_temp.csv"
                    filename = filename.replace(os.sep, '/')
                    data.to_csv(filename, sep=";", header=True, index=False)

                    # Create log files containing treated files level 0-1A-temp
                    if os.path.isfile(self.root + "/log/logfileDO_01A_temp.txt"):
                        write_option = 'a'
                    else:
                        with open(self.root + "/log/logfileDO_01A_temp.txt", "w") as f:
                            f.close()

                    if filename in open(self.root + "/log/logfileDO_01A_temp.txt").read():
                        next
                    else:
                        with open(self.root + "/log/logfileDO_01A_temp.txt", "a") as f:
                            f.write(filename + "\n")
                            f.close()

    def DO_maintenance_extract():
        filename = "../renku_ppmooring/ppmooring/scripts/PAR_maintenance.csv"
        maintenance = pd.read_csv(filename, sep=";")

        vec = np.copy(maintenance.debut)
        min_interval = min(abs(np.diff(vec)))  # smallest time interval

        minvec = min(vec)
        maxvec = max(vec)

        vecnew = np.arange(minvec, maxvec, min_interval)

        # Find index of maintenance points
        xy, xind, yind = np.intersect1d(vec, vecnew, return_indices=True)

        # y=np.ones(len(vecnew))
        # plt.plot(np.arange(len(vecnew)),y)
        # plt.plot(yind,y[yind],'.k')

        # Create vec_index
        vec_index = np.zeros(len(vecnew))
        vec_index[yind] = 1

        # index of the end of each interval
        end = []
        for i in range(0, (len(vec_index) - 1)):
            val = vec_index[i] - vec_index[i + 1]
            if val == 1:
                end = np.append(end, i)

        end = end.astype(int)
        end = vecnew[end]

        # index of the start of each interval
        debut = np.array([0])
        for i in range(0, (len(vec_index) - 1)):
            val = vec_index[i] - vec_index[i + 1]
            if val == -1:
                debut = np.append(debut, i)
        debut = debut.astype(int)
        debut = vecnew[debut]

        # Create dataframe
        maintenance = {"debut": debut, "end": end}
        maintenance = pd.DataFrame(maintenance)

        return maintenance

    #Quality flags in Level 0
    def quality_flags_level0(self, save=False):
        ## Perform quality assurance tests on data
        ## Code quality: 0 = nothing to report, 1 = more investigation

        # negative values = 0:
        ind = np.where(self.O2 < 0)[0]
        if len(ind) != 0:
            self.O2[ind] = 10000

        vec = np.copy(self.O2)
        vec_qual = np.copy(self.O2_qual)

        # Define variables
        variables = ["O2", "O2Sat"]

        # Remove maintenance data
        #maintenance = DO_maintenance_extract()

        indanom = []  # index anomaly
        for a in range(0, len(self.maintenance["debut"])):
            ind0 = np.where(self.time >= self.maintenance["debut"][a])[0]
            ind1 = np.where(self.time <= self.maintenance["end"][a])[0]
            ind = np.intersect1d(ind0, ind1)

            if len(ind) != 0:
                indanom = np.append(indanom, ind)
                indanom = np.unique(indanom, axis=0)

        # First filter:
        # remove maintenance dates:
        if len(indanom) != 0:
            vec_qual[indanom.astype(int)] = 1

        # nan values
        ind = np.where(np.isnan(vec) == True)[0]
        if len(ind) != 0:
            vec[ind] = 10000

        # Remove values > 20 mg/L
        ind = np.where(vec >= 20)[0]
        if len(ind) != 0:
            vec_qual[ind] = 1

        # first n values and last n values
        # first n values
        n = 30
        y_sub = np.copy(vec[0:n])
        y_qual = np.zeros(len(y_sub))
        y_qual = rm_outliers.rm_IQR(y_sub, y_qual, factor=3)
        if sum(y_qual)>0:
            vec_qual[0:n][y_qual] = 1  # update vec_qual

        # Last n values
        y_sub = np.copy(vec[-n:])
        y_qual = np.zeros(len(y_sub))
        y_qual = rm_outliers.rm_IQR(y_sub, y_qual, factor=3)
        if sum(y_qual)>0:
            vec_qual[-n:][y_qual] = 1  # update vec_qual

        # Remove threshold values error
        vdepth = self.depth[0]

        valid_range = {
            "50": [4, 20],
            "250": [4, 20],
            "500": [4, 20],
            "1000": [4, 20],
            "1500": [4, 20],
            "2000": [4, 15],
            "3000": [4, 15],
            "5000": [4, 15],
            "10000": [4, 15],
        }

        # Check if variable is within valid range

        indlow = np.where(vec < valid_range[str(vdepth)][0])
        indhigh = np.where(vec >= valid_range[str(vdepth)][1])
        vec_qual[indlow] = 1
        vec_qual[indhigh] = 1

        # Update vector O2_qual
        self.O2_qual = np.copy(vec_qual)

        # o2 anomaly = o2sat anomaly
        o2anomaly = np.where(self.O2_qual >= 1)[0]
        self.O2Sat_qual[o2anomaly] = 1

    #Quality flags in Level 1A
    def quality_flags_level1A(self, save=False):
        ## Perform quality assurance tests on data
        ## Code quality: 0 = nothing to report, 1 = more investigation

        # Variable O2, because O2_sat is calculated from O2
        vec = np.copy(self.O2)
        vec_qual = np.copy(self.O2_qual)
        id_point = np.array(range(0, len(vec)))

        # remove points precedent filter
        ind = np.where(vec_qual >= 1)[0]
        id_point = np.delete(id_point, ind)

        # Reinitialize
        vec = np.copy(self.O2)[id_point]
        vec_qual = np.zeros(len(vec))

        # Define valid range for each variable (min, max)
        border_val = pd.read_csv("../renku_ppmooring/ppmooring/scripts/border_value_lower.csv",sep=(';'))
        border_val[["min_depth", "min", "max"]] = border_val[["min_depth", "min", "max"]].apply(pd.to_numeric)
        border_val["min_depth"] = border_val["min_depth"] * 100

        # Depending on depth, value O2 can be relevant to be anormal
        vdepth = self.depth[0]
        vmonth = datetime.utcfromtimestamp(self.time[0]).month
        ind0 = np.where(border_val["month"] == vmonth)[0]
        ind1 = np.where(border_val["min_depth"] == vdepth)[0]
        if len(ind1) == 0:
            ind1 = np.where(border_val["min_depth"] == 0)[0]

        ind = np.intersect1d(ind0, ind1)

        vo2a = np.floor(border_val["min"][ind] * 0.95)
        vo2b = np.round(border_val["max"][ind] * 1.2)

        valid_range = {
            "O2": [vo2a, vo2b],
            "O2Sat": [0, 160]
        }
        # Check if variable is within valid range
        valid_values = (vec < float(valid_range["O2"][0])) | (vec > float(valid_range["O2"][1]))
        inf_values = np.isinf(vec)

        # Store to DataFrame
        vec_qual = (valid_values | inf_values).astype('uint8')

        # For O2_sat
        valid_values = (self.O2Sat < float(valid_range["O2Sat"][0])) | (self.O2 > float(valid_range["O2Sat"][1]))
        inf_values = np.isinf(self.O2Sat)
        self.O2Sat_qual = (valid_values | inf_values).astype('uint8')

        # remove points precedent filter
        ind = np.where(vec_qual >= 1)[0]
        id_point = np.delete(id_point, ind)

        # Reinitialize
        vec = np.copy(self.O2)[id_point]
        vec_qual = np.zeros(len(vec))

        # Second filer: variation rate between two successive values.
        # rm variation rate
        if len(vec) != 0:
            y_qual = rm_outliers.rm_variation_rate(vec, vec_qual)
        else:
            y_qual = np.zeros(len(vec))

        if sum(y_qual)>0:
            id_point = np.delete(id_point, y_qual)
            y_qual =[]

            # Reinitialize
            vec = np.copy(self.O2)[id_point]
            vec_qual = np.zeros(len(vec))

        if len(vec)<10:
            # Final: all values superior than 1 in vec_qual become 1
            self.O2_qual = np.zeros(len(self.O2_qual)) + 1
            self.O2_qual[id_point] = 0

            # o2 anomaly = o2sat anomaly
            o2anomaly = np.where(self.O2_qual >= 1)[0]
            self.O2Sat_qual[o2anomaly] = 1
            return

        # #Third filter:
        # Configure
        from scipy.signal import savgol_filter
        import pwlf

        n_changepoints = 1
        smoothing_window_size = 31
        smoothing_order = 3

        # Prepare (smoothed) data
        x = np.array(range(0, len(vec)))
        vec_smooth = savgol_filter(vec, smoothing_window_size, smoothing_order)

        # Fit PLR
        plr_model = pwlf.PiecewiseLinFit(x, vec_smooth)
        breaks = plr_model.fit(n_changepoints + 1)
        x_hat = np.linspace(x.min(), x.max(), 100)
        y_hat = plr_model.predict(x_hat)

        # reconstruction
        if len(breaks) < 3:
            ibreak = breaks[1].astype(int)
            mean_before = np.median(vec_smooth[0:ibreak])
            mean_after = np.median(vec_smooth[ibreak:len(vec_smooth)])

            gap = mean_after - mean_before

            vecnew = np.copy(vec)
            if len(vec_smooth[0:ibreak]) < len(
                    vec_smooth[ibreak:len(vec_smooth)]):  # length of first segment < second segment
                if gap > 0:
                    vecnew[0:ibreak] = vecnew[0:ibreak] + gap
                    i_reconstruct = id_point[0:ibreak]
                else:
                    vecnew[0:ibreak] = vecnew[0:ibreak] - gap
                    i_reconstruct = id_point[0:ibreak]
            else:
                if gap > 0:
                    vecnew[ibreak:len(vecnew)] = vecnew[ibreak:len(vecnew)] - gap
                    i_reconstruct = id_point[ibreak:len(vecnew)]
                else:
                    vecnew[ibreak:len(vecnew)] = vecnew[ibreak:len(vecnew)] + gap
                    i_reconstruct = id_point[ibreak:len(vecnew)]

            vec = np.copy(vecnew)

            # reconstruct
            if len(i_reconstruct) != 0:
                self.O2[i_reconstruct] = np.copy(vec[i_reconstruct])

        if len(breaks) >= 3:
            y_qual = rm_outliers.rm_kmeans_threshold(vec, ncluster=2,threshold=0.2)

        if sum(y_qual)>0:
            id_point = np.delete(id_point, y_qual)
            #idx = []
            # Reinitialize
            vec = np.copy(self.O2)[id_point]
            y_qual = []
            vec_qual = np.zeros(len(vec))

        #Next _filter
        xval = max(vec) / np.quantile(vec, 0.95)
        if xval > 1.1:
            y_qual = rm_outliers.rm_kmeans_threshold(vec, ncluster=2, threshold=0.1)
            if sum(y_qual)>0:
                id_point = np.delete(id_point, y_qual)
                #idx = []
                # Reinitialize
                vec = np.copy(self.O2)[id_point]
                vec_qual = np.zeros(len(vec))
                y_qual = []

        # Next _filter
        xval = max(vec) / np.quantile(vec, 0.95)
        if xval > 1.05:
            dayperiode = (24 * 60 * 60 / (self.time[1] - self.time[0])).astype(int)
            if len(vec) >= 3 * dayperiode:
                y_qual = rm_outliers.rm_convolution(vec, vec_qual, window_len=dayperiode, window_type='blackman',
                                             threshold=np.quantile(vec, 0.99), n_sigma=2.5)
            else:
                y_qual = rm_outliers.rm_convolution(vec, vec_qual, window_len=round(len(vec) / 2), window_type='blackman',
                                             threshold=np.quantile(vec, 0.99), n_sigma=2.5)
        if sum(y_qual)>0:
            id_point = np.delete(id_point, y_qual)
            y_qual =[]
            # Reinitialize
            vec = np.copy(self.O2)[id_point]
            vec_qual = np.zeros(len(vec))

        # Next _filter
        xval = max(vec) - np.quantile(vec, 0.95)
        if xval > 0.5:
            idx = np.where(vec >= np.quantile(vec, 0.999))[0]

            if len(idx) != 0:
                id_point = np.delete(id_point, idx)
                y_qual = []
                # Reinitialize
                vec = np.copy(self.O2)[id_point]
                vec_qual = np.zeros(len(vec))

        xval = max(vec) - np.quantile(vec, 0.95)
        if xval > 0.5:
            y_qual = rm_outliers.rm_variation_rate(vec, vec_qual)

        if sum(y_qual)>0:
            id_point = np.delete(id_point, y_qual)

            # Reinitialize
            vec = np.copy(self.O2)[id_point]
            vec_qual = np.zeros(len(vec))
            y_qual =[]

        #Last filter
        if len(vec) >= 5:
            y_qual = rm_outliers.rm_max(vec, vec_qual, factor1=30, semiwindow=round(len(vec) / 2))
        else:
            idx = []

        if sum(y_qual)>0:
            id_point = np.delete(id_point, y_qual)
            y_qual = []
            vec_qual = np.zeros(len(vec))

        # Final: all values superior than 1 in vec_qual become 1
        self.O2_qual = np.zeros(len(self.O2_qual)) + 1
        self.O2_qual[id_point] = 0

        # o2 anomaly = o2sat anomaly
        o2anomaly = np.where(self.O2_qual >= 1)[0]
        self.O2Sat_qual[o2anomaly] = 1

class PAR_data(ppmooring):
    ## Class for the PAR data

    def __init__(self):
        ## Set attributes for this class
        self.root = ppmooring.ppmooringroot
        # Output folders
        self.outfolders = {
            "RBR_PAR_50cm": "/LexplorePPMooringPAR_50cm",
            "RBR_PAR_250cm": "/LexplorePPMooringPAR_250cm",
            "RBR_PAR_500cm": "/LexplorePPMooringPAR_500cm",
            "RBR_PAR_1000cm": "/LexplorePPMooringPAR_1000cm",
            "RBR_PAR_2000cm": "/LexplorePPMooringPAR_2000cm",
            "RBR_PAR_3000cm": "/LexplorePPMooringPAR_3000cm",
            "50cm": "/LexplorePPMooringPAR_50cm",
            "250cm": "/LexplorePPMooringPAR_250cm",
            "500cm": "/LexplorePPMooringPAR_500cm",
            "1000cm": "/LexplorePPMooringPAR_1000cm",
            "2000cm": "/LexplorePPMooringPAR_2000cm",
            "3000cm": "/LexplorePPMooringPAR_3000cm"
        }
        self.Level2_outfolder = "/LexplorePPMooringPAR"

        # Matrix depth - folder
        self.depthdict = {
            "RBR_PAR_50cm":50,
            "RBR_PAR_250cm":250,
            "RBR_PAR_500cm":500,
            "RBR_PAR_1000cm":1000,
            "RBR_PAR_2000cm":2000,
            "RBR_PAR_3000cm":3000
        }

        # General attributes
        self.att = {
            "institution": "EPFL",
            "source": "Mooring DO",
            "references": "LéXPLORE commun instruments Viet Tran Khac <viet.tran-khac@inrae.fr>",
            "history": "See history on Renku",
            "conventions": "CF 1.7",
            "comment": "PAR (photosynthetically active radiation) from Mooring M1 on Lexplore Platform in Lake Geneva",
            "title": "Lexplore PPMooring PAR"
        }

        # Dimension attributes
        self.dim_dict = {
            'time': {'dim_name':'time', 'dim_size':None}
        }

        # Coordinate variable attributes
        self.dim_var_dict = {
            'time': {'var_name':'time', 'dim':('time',), 'unit':'seconds since 1970-01-01 00:00:00',
                     'longname':'time'}
        }

        # Level 1A - Variable attributes
        self.L1A_dict = {
            'depth': {'var_name':'depth', 'unit':'cm', 'longname':'depth'},
            'PAR': {'var_name':'PAR', 'dim':('time',), 'unit':'μmol/m²/s', 'longname':'PAR'},
            'PAR_qual': {'var_name': 'PAR_qual', 'dim': ('time',),
                         'unit': '0 = nothing to report, 1 = more investigation',
                         'longname': 'PAR quality flag'}

        }

        # Level 1B - Quality flag attributes
        self.L1B_dict = {
            'PAR_qual': {'var_name':'PAR_qual', 'dim':('time',), 'unit':'0 = nothing to report, 1 = more investigation',
                         'longname':'PAR quality flag'}
        }

        # Level 2 - Variable attributes
        self.L2_dict = {
            'PAR': {'var_name':'PAR', 'dim':('time',), 'unit':'μmol/m²/s', 'longname':'PAR'}
        }

    def get_files_to_read_Level0(self, input_dir):
        ## Get a list of all files to read

        # List all subfolders
        folders = os.listdir(input_dir)
        folders = fnmatch.filter(folders, 'RBR_*')
        folders = [input_dir + "/" + s + "/conv" for s in folders]

        self.filestoread = []
        # Loop through all subfolders
        for idxfolder, folder in enumerate(folders):
             files = os.listdir(folder)
             files_subfolder_txt = [folder + "/" + s for s in fnmatch.filter(files, "*.txt")]
             files_subfolder_xlsx = [folder + "/" + s for s in fnmatch.filter(files, "*.xlsx")]
             files_subfolder = files_subfolder_txt + files_subfolder_xlsx
             self.filestoread = np.append(self.filestoread, files_subfolder)

    def get_files_to_read_Leveltemp(self, input_dir):
        ## Get a list of all files to read at level n

        # List all subfolders
        folders = os.listdir(input_dir)
        folders = fnmatch.filter(folders, 'LexplorePPMooringPAR*')
        folders = [input_dir + "/" + s + "/temp" for s in folders]

        self.filestoread = []
        # Loop through all subfolders
        for idxfolder, folder in enumerate(folders):
            files = os.listdir(folder)
            files_subfolder = [folder + "/" + s for s in fnmatch.filter(files, "*.csv")]
            self.filestoread = np.append(self.filestoread, files_subfolder)

    def read_data(self, infile):
        #Since update of software RBR version 6.13, data is converted into .xlsx, no longer .txt

        #Read txt files
        if (infile.endswith('txt')):
            ## Read PAR data

            # Store infile variable
            self.infile = infile

            # Get input folder
            self.folder = self.infile.rsplit("/")[-3]

            # Import metadata and get indicies of metadata
            with open(self.infile, encoding="utf8", errors='ignore') as f:
                for idx, line in enumerate(f):
                    if 'Serial' in line:
                        idx_serial = idx
                    if 'Channel[1].calibration' in line:
                        idx_calibr = idx
                    if 'NumberOfSamples' in line:
                        idx_meta = idx + 1
                        break
                f.seek(0)
                lines = f.readlines()
                meta_data = lines[0:idx_meta]

            # Extract meta-data
            numsensor = meta_data[idx_serial].rsplit("=")[1].rstrip("\n")
            numcalibration = meta_data[idx_calibr].rsplit("=")[1].rstrip("\n")

            # Read Data
            df = pd.read_csv(self.infile, delim_whitespace=True, header=None, skiprows=idx_meta + 2,encoding="unicode_escape")
            df.columns = ["Date", "Time", "PAR"]

            # Replace date and time by timestamp
            datehour = df["Date"] + " " + df["Time"]
            datehour = datehour.apply(lambda x: datetime.strptime(x, '%Y-%m-%d %H:%M:%S.%f'))
            df["Timestamp"] = [datetime.timestamp(s) for s in datehour]

            df = df[["Timestamp", "PAR"]]

            # Sort by date and reset index
            df = df.sort_values(by=['Timestamp'])
            self.df = df.reset_index(inplace=False, drop=True)

            # Sensor data, depth and point
            self.meta_data = {
                "sensor": numsensor,
                "calibration": numcalibration,
                "instrument": "rbr",
                "depth": self.depthdict[self.folder],
                "point": self.folder
            }

            # Store variables
            self.time = df["Timestamp"]  # need to be checked
            self.depth = self.meta_data["depth"]
            self.PAR = df["PAR"]

            # Create vector quality for data quality flags at Level 1A
            self.PAR_qual = np.zeros(len(self.PAR))

        #Read xlsx files
        if (infile.endswith('xlsx')):
            ## Read PAR data

            # Store infile variable
            self.infile = infile

            # Get input folder
            self.folder = self.infile.rsplit("/")[-3]

            # Import metadata
            metadata = pd.read_excel(open(infile, 'rb'),sheet_name = 'Metadata')

            #Get metadata
            metadata_nrow = metadata.shape[0]
            metadata_ncol = metadata.shape[1]

            idx = np.where(metadata.iloc[:,0].str.contains("Serial")==True)[0][0] +1
            numsensor = metadata.iloc[idx,0]

            idx = np.where(metadata.iloc[:,0].str.contains("Atmospheric pressure")==True)[0][0] +1
            atmospheric_pressure = metadata.iloc[idx,0]
            density = metadata.iloc[idx,1]
            Altitude = metadata.iloc[idx,2]
            meta_Temperature = metadata.iloc[idx,3]
            pressure = metadata.iloc[idx,4]
            meta_conductivity = metadata.iloc[idx,5]
            salinity = metadata.iloc[idx,6]
            conductivity_specific = metadata.iloc[idx,7]


            # Read Data
            df = pd.read_excel(open(infile, 'rb'),sheet_name = 'Data',skiprows=1)
            df.columns = ["Time", "PAR"]

            # Replace date and time by timestamp
            df["Timestamp"] = [s.timestamp() for s in df["Time"]]

            df = df[["Timestamp", "PAR"]]

            # Sort by date and reset index
            df = df.sort_values(by=['Timestamp'])
            self.df = df.reset_index(inplace=False, drop=True)

            # Sensor data, depth and point
            self.meta_data = {
                "sensor": numsensor,
                "atmospheric_pressure" : atmospheric_pressure,
                "density": density,
                "Altitude": Altitude,
                "meta_Temperature": meta_Temperature,
                "pressure": pressure,
                "meta_conductivity": meta_conductivity,
                "salinity": salinity,
                "conductivity_specific": conductivity_specific,
                "instrument": "rbr",
                "depth": self.depthdict[self.folder],
                "point": self.folder
            }

            # Store variables
            self.time = df["Timestamp"]  # need to be checked
            self.depth = self.meta_data["depth"]
            self.PAR = df["PAR"]

            # Create vector quality for data quality flags at Level 1A
            self.PAR_qual = np.zeros(len(self.PAR))

    def read_data_Leveln(self, infile):
        # Read data
        df = pd.read_csv(infile, sep=";")  # delim_whitespace=True, header=None, skiprows=idx_meta + 2)
        df = df.reset_index(drop=True)

        # Replace date and time by timestamp
        df["Timestamp"] = df["time"]

        df = df[["Timestamp", "PAR", "PAR_qual","depth"]]

        # Sort by date and reset index
        df = df.sort_values(by=['Timestamp'])
        self.df = df.reset_index(inplace=False, drop=True)

        # Store variables
        # Get input folder
        self.name1B = infile
        self.folder = str(df["depth"][0]) + "cm"

        self.time = df["Timestamp"]  # need to be checked
        self.depth = df["depth"]
        self.PAR = df["PAR"]
        self.PAR_qual = df["PAR_qual"]

    def parsetimedata(self, atlevel, freqday):
        root = self.root
        input_dir = os.path.join(root, "data", atlevel)
        input_dirs = os.listdir(input_dir)
        input_dirs = fnmatch.filter(input_dirs, 'LexplorePPMooringPAR_*')

        for folder in input_dirs:
            depthfolder = os.path.join(input_dir, folder)
            files = os.listdir(depthfolder)
            files = fnmatch.filter(files,'*.csv')
            ##files = files[0:(len(files) - 1)]  # remove temp folder

            dfall = pd.DataFrame()
            for file in files:
                filepath = os.path.join(depthfolder, file)
                df = pd.read_csv(filepath)
                #dfall = dfall_append(df) #pandas <2.0
                dfall = pd.concat([dfall, df], ignore_index=True)

            # Split into 15 days files
            dfall['time1'] = dfall.time.apply(lambda x: datetime.fromtimestamp(x))

            # Split df_final into 15 day dataframes and expor temporarily into csv
            # freqday = '15D'
            days = [g for n, g in dfall.groupby(pd.Grouper(key="time1", freq=freqday))]

            for day in range(0, len(days)):
                if len(days[day]) == 0:
                    next
                else:
                    data = days[day]
                    data = data.sort_values(by=['time'])
                    data = data.reset_index(inplace=False, drop=True)

                    outfolder = os.path.join(depthfolder, "temp")
                    if not os.path.exists(outfolder):
                        os.makedirs(outfolder)

                    prefix_start = datetime.fromtimestamp(data["time"].iloc[0]).strftime('%Y%m%d')
                    prefix_end = datetime.fromtimestamp(data["time"].iloc[-1]).strftime('%Y%m%d')

                    datename = prefix_start + "_" + prefix_end

                    filename = outfolder + "/" + datename + "_" + folder + "Level1A_temp.csv"
                    filename = filename.replace(os.sep, '/')
                    data.to_csv(filename, sep=";", header=True, index=False)

                    # Create log files containing treated files level 0-1A-temp
                    if os.path.isfile(self.root + "/log/logfile_01A_temp.txt"):
                        write_option = 'a'
                    else:
                        with open(self.root + "/log/logfile_01A_temp.txt", "w") as f:
                            f.close()

                    if filename in open(self.root + "/log/logfile_01A_temp.txt").read():
                        next
                    else:
                        with open(self.root + "/log/logfile_01A_temp.txt", "a") as f:
                            f.write(filename + "\n")
                            f.close()

        # Quality flags in Level 0

    def quality_flags_level0(self, save=False):
        ## Perform quality assurance tests on data
        ## Code quality: 0 = nothing to report, 1 = more investigation

        # First filter: at the beginning and the end of the vector

        # negative values = 0:
        ind = np.where(self.PAR < 0)[0]
        if len(ind) != 0:
            self.PAR[ind] = 0

        vec = np.copy(self.PAR)
        vec_qual = np.copy(self.PAR_qual)


        maintenance = pd.read_csv(self.root + "/scripts/PAR_maintenance.csv", sep=";")
        indanom = []  # index anomaly
        for a in range(0, len(maintenance["debut"])):
            # ind = np.where(self.time >= maintenance["debut"][a])[0]
            # ind1 = np.where(self.time <= maintenance["end"][a])[0]
            # ind2 = list(set(ind) & set(ind1))
            ind = np.where(self.time == maintenance["debut"][a])[0]

            if len(ind) != 0:
                indanom = np.append(indanom, ind)
                indanom = np.unique(indanom,axis=0)

        # First filter:
        # remove maintenance dates:
        if len(indanom) != 0:
            vec_qual[indanom.astype(int)] = 1

        # nan values
        ind = np.where(np.isnan(vec)==True)[0]
        if len(ind) != 0:
            vec[ind] = 10000

        #Remove values > 10000 μmol s-1 m-2
        ind = np.where(vec>=10000)[0]
        if len(ind) !=0:
            vec_qual[ind] = 1

        # first n values and last n values
        # first n values
        n = 30
        y_sub = np.copy(vec[0:n])
        y_qual = np.zeros(len(y_sub))
        y_qual = rm_outliers.rm_IQR(y_sub, y_qual, factor=3)

        if sum(y_qual)>0:
            vec_qual[0:n][y_qual] = 1  # update vec_qual

        # Last n values
        y_sub = np.copy(vec[-n:])
        y_qual = np.zeros(len(y_sub))

        y_qual = rm_outliers.rm_IQR(y_sub, y_qual, factor=3)
        if sum(y_qual)>0:
            vec_qual[-n:][y_qual] = 1  # update vec_qual

        #Remove threshold values error
        vdepth = self.depth

        valid_range = {
            "50": [0, 2000],
            "250": [0, 2000],
            "500": [0, 1500],
            "1000": [0, 1000],
            "2000": [0, 30],
            "3000": [0, 10]}

        # Check if variable is within valid range
        if isinstance(vdepth,int)==True:
            indlow = np.where(vec<valid_range[str(vdepth)][0])
            indhigh = np.where(vec>=valid_range[str(vdepth)][1])
        else:
            indlow = np.where(vec < valid_range[str(vdepth[0])][0])
            indhigh = np.where(vec >= valid_range[str(vdepth[0])][1])
        vec_qual[indlow] = 1
        vec_qual[indhigh] = 1

        # Update vector PAR_qual
        self.PAR_qual = np.copy(vec_qual)

    # Quality flags in Level 1A
    def quality_flags_level1A(self, save=False):
        ## Perform quality assurance tests on data
        ## Code quality: 0 = nothing to report, 1 = more investigation

        vec = np.copy(self.PAR)
        vec_qual = np.copy(self.PAR_qual)
        id_point = np.array(range(0, len(vec)))

        # remove points precedent filter
        ind = np.where(vec_qual >= 1)[0]
        id_point = np.delete(id_point, ind)

        #Reinitialize
        vec = np.copy(self.PAR)[id_point]
        vec_qual = np.zeros(len(vec))

        #if vector is empty = file does not contain data
        if len(vec)!=0:
            #1st filter: kmean cluster
            y_qual = rm_outliers.rm_kmeans(vec,ncluster=2)

            if sum(y_qual)>0:
                id_point = np.delete(id_point, y_qual)

        # Reinitialize
        vec = np.copy(self.PAR)[id_point]
        vec_qual = np.zeros(len(vec))
        y_qual =[]

        #2nd filter: rm rm_convolution
        dayperiode = (24 * 60 * 60 / (self.time[1] - self.time[0])).astype(int)

        if len(vec) >= 3 * dayperiode:
            y_qual = rm_outliers.rm_convolution(vec, vec_qual, window_len=dayperiode, window_type='blackman',
                                         threshold=np.quantile(vec, 0.999), n_sigma=2.5)
        else:
            y_qual = rm_outliers.rm_convolution(vec, vec_qual, window_len=round(dayperiode / 2), window_type='blackman',
                                         threshold=np.quantile(vec, 0.999), n_sigma=2.5)

        if sum(y_qual)>0:
            id_point = np.delete(id_point, y_qual)

        # Reinitialize
        vec = np.copy(self.PAR)[id_point]
        vec_qual = np.zeros(len(vec))
        y_qual = []

        #3rd filter: rm variation rate
        y_qual = rm_outliers.rm_variation_rate(vec,vec_qual)

        if sum(y_qual)>0:
            id_point = np.delete(id_point, y_qual)

        #Update vector PAR_qual
        self.PAR_qual = np.zeros(len(self.PAR_qual))+1
        self.PAR_qual[id_point] = 0

class Temperature_data(ppmooring):
    ## Class for the Temperature data

    def __init__(self):
        ## Set attributes for this class
        self.root = ppmooring.ppmooringroot
        # Output folders
        self.outfolders = {
            "T_1": "/LexplorePPMooringTemperature_50cm",
            "T_2": "/LexplorePPMooringTemperature_250cm",
            "T_3": "/LexplorePPMooringTemperature_500cm",
            "T_4": "/LexplorePPMooringTemperature_750cm",
            "T_5": "/LexplorePPMooringTemperature_1000cm",
            "T_6": "/LexplorePPMooringTemperature_1250cm",
            "T_7": "/LexplorePPMooringTemperature_1500cm",
            "T_8": "/LexplorePPMooringTemperature_1750cm",
            "T_9": "/LexplorePPMooringTemperature_2000cm",
            "T_10": "/LexplorePPMooringTemperature_2250cm",
            "T_11": "/LexplorePPMooringTemperature_2500cm",
            "T_12": "/LexplorePPMooringTemperature_2750cm",
            "T_13": "/LexplorePPMooringTemperature_3000cm",
            "50cm": "/LexplorePPMooringTemperature_50cm",
            "250cm": "/LexplorePPMooringTemperature_250cm",
            "500cm": "/LexplorePPMooringTemperature_500cm",
            "750cm": "/LexplorePPMooringTemperature_750cm",
            "1000cm": "/LexplorePPMooringTemperature_1000cm",
            "1250cm": "/LexplorePPMooringTemperature_1250cm",
            "1500cm": "/LexplorePPMooringTemperature_1500cm",
            "1750cm": "/LexplorePPMooringTemperature_1750cm",
            "2000cm": "/LexplorePPMooringTemperature_2000cm",
            "2250cm": "/LexplorePPMooringTemperature_2250cm",
            "2500cm": "/LexplorePPMooringTemperature_2500cm",
            "2750cm": "/LexplorePPMooringTemperature_2750cm",
            "3000cm": "/LexplorePPMooringTemperature_3000cm"
        }
        self.Level2_outfolder = "/LexplorePPMooringTemperature"

        # Matrix depth - folder
        self.depthdict = {
            "T_1":50,
            "T_2":250,
            "T_3":500,
            "T_4":750,
            "T_5":1000,
            "T_6":1250,
            "T_7":1500,
            "T_8":1750,
            "T_9":2000,
            "T_10":2250,
            "T_11":2500,
            "T_12":2750,
            "T_13":3000
        }

        # General attributes
        self.att = {
            "institution": "EPFL",
            "source": "PPMooring DO",
            "references": "LéXPLORE commun instruments Viet Tran Khac <viet.tran-khac@inrae.fr>",
            "history": "See history on Renku",
            "conventions": "CF 1.7",
            "comment": "Water temperatures from PPMooring M2 on Lexplore Platform in Lake Geneva",
            "title": "Lexplore PPMooring Temperature"
        }

        # Dimension attributes
        self.dim_dict = {
            'time': {'dim_name':'time', 'dim_size':None}
        }

        # Coordinate variable attributes
        self.dim_var_dict = {
            'time': {'var_name':'time', 'dim':('time',), 'unit':'seconds since 1970-01-01 00:00:00', 'longname':'time'}
        }

        # Level 1A - Variable attributes
        self.L1A_dict = {
            'depth': {'var_name':'depth', 'unit':'cm', 'longname':'depth'},
            'Temp': {'var_name':'Temp', 'dim':('time',), 'unit':'degC', 'longname':'temperature'},
            'Temp_qual': {'var_name': 'Temp_qual', 'dim': ('time',),
                          'unit': '0 = nothing to report, 1 = more investigation',
                          'longname': 'temperature quality flag'}
        }

        # Level 1B - Quality flag attributes
        self.L1B_dict = {
            'Temp_qual': {'var_name':'Temp_qual', 'dim':('time',), 'unit':'0 = nothing to report, 1 = more investigation', 'longname':'temperature quality flag'}
        }

        # Level 2 - Variable attributes
        self.L2_dict = {
            'Temp': {'var_name':'Temp', 'dim':('time',), 'unit':'degC', 'longname':'temperature'}
        }

    def get_files_to_read_Level0(self, input_dir):
        ## Get a list of all files to read

        # List all subfolders
        folders = os.listdir(input_dir)
        folders = fnmatch.filter(folders, 'T_*')
        folders = [input_dir + "/" + s + "/conv" for s in folders]

        self.filestoread = []
        # Loop through all subfolders
        for idxfolder, folder in enumerate(folders):
            files = os.listdir(folder)
            self.filestoread = np.append(self.filestoread, [folder + "/" + s for s in fnmatch.filter(files, "*.txt")])
            self.filestoread = np.append(self.filestoread, [folder + "/" + s for s in fnmatch.filter(files, "*.csv")])


    def get_files_to_read_Leveltemp(self, input_dir):
        ## Get a list of all files to read

        # List all subfolders
        folders = os.listdir(input_dir)
        folders = fnmatch.filter(folders, '*Temperature*')
        folders = [input_dir + "/" + s + "/temp" for s in folders]

        self.filestoread = []
        # Loop through all subfolders
        for idxfolder, folder in enumerate(folders):
            files = os.listdir(folder)
            files_subfolder = [folder + "/" + s for s in fnmatch.filter(files, "*.csv")]
            self.filestoread = np.append(self.filestoread, files_subfolder)

    def read_data_txt(self, infile):
        ## Read .txt data

        # Store infile variable
        self.infile = infile

        # Get input folder
        self.folder = self.infile.rsplit("/")[-3]

        # Import metadata and get indicies of metadata
        with open(self.infile, encoding="utf8", errors='ignore') as f:
            for idx, line in enumerate(f):
                if 'Serial' in line:
                    idx_serial = idx
                if 'Channel[1].calibration' in line:
                    idx_calibr = idx
                if 'NumberOfSamples' in line:
                    idx_meta = idx + 1
                    break
            f.seek(0)
            lines = f.readlines()
            meta_data = lines[0:idx_meta]

        # Extract meta-data
        numsensor = meta_data[idx_serial].rsplit("=")[1].rstrip("\n")
        numcalibration = meta_data[idx_calibr].rsplit("=")[1].rstrip("\n")

        # Read data
        df = pd.read_csv(self.infile, delim_whitespace=True, header=None, skiprows=idx_meta + 2)
        df.columns = ["Date", "Time", "Temp"]

        # Replace date and time by timestamp
        datehour = [df["Date"] + " " + df["Time"]]
        datehour = datehour[0].apply(lambda x: datetime.strptime(x, '%d-%b-%Y %H:%M:%S.%f'))
        df["Timestamp"] = [datetime.timestamp(s) for s in datehour]

        df = df[["Timestamp", "Temp"]]

        # Sort by date and reset index
        df = df.sort_values(by=['Timestamp'])
        self.df = df.reset_index(inplace=False, drop=True)

        # Sensor data, depth and point
        self.meta_data = {
            "sensor": numsensor,
            "calibration": numcalibration,
            "instrument": "rbr",
            "depth": self.depthdict[self.folder],
            "point": self.folder
        }

        # Store variables
        self.time = df["Timestamp"]  # need to be checked
        self.depth = self.meta_data["depth"]
        self.Temp = df["Temp"]

        # Create vector quality for data quality flags at Level 1A
        self.Temp_qual = np.zeros(len(self.Temp))

    def read_data_csv(self, infile):
        ## Read .txt data

        # Store infile variable
        self.infile = infile

        # Get input folder
        self.folder = self.infile.rsplit("/")[-3]

        # Import metadata and get indicies of metadata
        with open(self.infile, encoding="utf8", errors='ignore') as f:
            for idx, line in enumerate(f):
                if 'Source Device' in line:
                    idx_srcdev = idx
                    break
            f.seek(0)
            lines = f.readlines()
            meta_data = lines[0:8]

        # Extract meta-data
        numsensor = meta_data[idx_srcdev].rsplit("-")[-1].rstrip("\n")
        numcalibration = -99

        # Read data
        df = pd.read_csv(self.infile, sep=",", header=None, skiprows=8)

        # Some files has ADC column, others do not have
        if len(df.columns) == 4:
            df = df.drop(columns=3)
        df.columns = ["Date", "Time", "Temp"]

        # Replace date and time by timestamp
        datehour = [df["Date"] + " " + df["Time"]]
        datehour = datehour[0].apply(lambda x: datetime.strptime(x, '%Y-%m-%d %H:%M:%S'))
        df["Timestamp"] = [datetime.timestamp(s) for s in datehour]

        df = df[["Timestamp", "Temp"]]

        # Sort by date and reset index
        df = df.sort_values(by=['Timestamp'])
        self.df = df.reset_index(inplace=False, drop=True)

        # Sensor data, depth and point
        self.meta_data = {
            "sensor": numsensor,
            "calibration": numcalibration,
            "instrument": "minidot",
            "depth": self.depthdict[self.folder],
            "point": self.folder
        }

        # Store variables
        self.time = df["Timestamp"]  # need to be checked
        self.depth = self.meta_data["depth"]
        self.Temp = df["Temp"]

        # Create vector quality for data quality flags at Level 1A
        self.Temp_qual = np.zeros(len(self.Temp))

    def read_data_Leveln(self, infile):
        # Read data
        df = pd.read_csv(infile, sep=";")  # delim_whitespace=True, header=None, skiprows=idx_meta + 2)
        df = df.reset_index(drop=True)

        # Replace date and time by timestamp
        df["Timestamp"] = df["time"]

        df = df[["Timestamp", "Temp", "Temp_qual", "depth"]]

        # Sort by date and reset index
        df = df.sort_values(by=['Timestamp'])
        self.df = df.reset_index(inplace=False, drop=True)

        # Store variables
        # Get input folder
        self.name1B = infile
        self.folder = str(df["depth"][0]) + "cm"

        self.time = df["Timestamp"]  # need to be checked
        self.depth = df["depth"]
        self.Temp = df["Temp"]
        self.Temp_qual = df["Temp_qual"]

    def parsetimedata(self, atlevel, freqday):
        root = self.root
        input_dir = os.path.join(root, "data", atlevel)
        input_dirs = os.listdir(input_dir)
        input_dirs = fnmatch.filter(input_dirs, '*Temperature*')

        for folder in input_dirs:
            depthfolder = os.path.join(input_dir, folder)
            files = os.listdir(depthfolder)
            files = fnmatch.filter(files, '*.csv')
            files = files[0:(len(files) - 1)]  # remove temp folder

            dfall = pd.DataFrame()
            for file in files:
                filepath = os.path.join(depthfolder, file)
                df = pd.read_csv(filepath)
                # dfall = dfall_append(df) #pandas <2.0
                dfall = pd.concat([dfall, df], ignore_index=True)

            # Split into 15 days files
            dfall['time1'] = dfall.time.apply(lambda x: datetime.fromtimestamp(x))

            # Split df_final into 15 day dataframes and expor temporarily into csv
            freqday = '15D'
            days = [g for n, g in dfall.groupby(pd.Grouper(key="time1", freq=freqday))]

            for day in range(0, len(days)):
                if len(days[day]) == 0:
                    next
                else:
                    data = days[day]
                    data = data.sort_values(by=['time'])
                    data = data.reset_index(inplace=False, drop=True)

                    outfolder = os.path.join(depthfolder, "temp")
                    if not os.path.exists(outfolder):
                        os.makedirs(outfolder)

                    prefix_start = datetime.fromtimestamp(data["time"].iloc[0]).strftime('%Y%m%d')
                    prefix_end = datetime.fromtimestamp(data["time"].iloc[-1]).strftime('%Y%m%d')

                    datename = prefix_start + "_" + prefix_end

                    filename = outfolder + "/" + datename + "_" + folder + "Level1A_temp.csv"
                    filename = filename.replace(os.sep, '/')
                    data.to_csv(filename, sep=";", header=True, index=False)

                    # Create log files containing treated files level 0-1A-temp
                    if os.path.isfile(self.root + "/log/logfile_01A_temp.txt"):
                        write_option = 'a'
                    else:
                        with open(self.root + "/log/logfile_01A_temp.txt", "w") as f:
                            f.close()

                    if filename in open(self.root + "/log/logfile_01A_temp.txt").read():
                        next
                    else:
                        with open(self.root + "/log/logfile_01A_temp.txt", "a") as f:
                            f.write(filename + "\n")
                            f.close()

        # Quality flags in Level 0
        def quality_flags_level0(self):
            ## Perform quality assurance tests on data
            ## Code quality: 0 = nothing to report, 1 = more investigation

            # Define variables
            variables = ["Temp"]

            # Define valid range for each variable (min, max)
            valid_range = {
                "temp":[0, 40]
            }

            # Loop over variables
            for var in variables:

                # Check if variable is within valid range
                if len(valid_range[var]) == 2:
                    valid_values = (getattr(self,var) < valid_range[var][0]) \
                        | (getattr(self,var) > valid_range[var][1])
                else:
                    valid_values = getattr(self,var) < valid_range[var][0]

                # Check if variable is a number
                datatype_values = np.full(getattr(self,var).shape, False)
                for idx, y in np.ndenumerate(getattr(self,var)):
                    try:
                        float(y)
                    except:
                        datatype_values[idx] = True

                # Check for NaN and Inf vlaues
                nan_values = np.isnan(getattr(self,var))
                inf_values = np.isinf(getattr(self,var))

                # Store to DataFrame
                setattr(self, var+'_qual', (valid_values | nan_values | inf_values | datatype_values).astype('uint8'))


        def quality_flags_level1A(self):
            ## Perform quality assurance tests on data
            ## Code quality: 0 = nothing to report, 1 = more investigation
            # All maintenance date = code quality = 1
            maintenance = pd.read_csv(self.root + "/scripts/PAR_maintenance.csv", sep=";")
            indanom = []  # index anomaly
            for a in range(0, len(maintenance["debut"])):
                ind = np.where(self.time >= maintenance["debut"][a])[0]
                ind1 = np.where(self.time <= maintenance["end"][a])[0]
                ind2 = list(set(ind) & set(ind1))

                if len(ind2) != 0:
                    indanom = np.append(indanom, ind2)

            # First filter:
            vec = np.copy(self.Temp)
            vec_qual = np.copy(self.Temp_qual)
            id_point = np.array(range(0, len(vec)))

            #negative values = 0:
            ind = np.where(vec<0)[0]
            self.Temp[ind] = 0
            vec = np.copy(self.Temp)

            #remove maintenance dates:
            if len(indanom)!=0:
                vec_qual[indanom.astype(int)] = 1

            # remove points precedent filter
            ind = np.where(vec_qual >= 1)[0]
            id_point = np.delete(id_point, ind)

            # reference is y_sub0 from now on, id_point is always updated
            vec = np.copy(self.Temp)[id_point]
            y_qual0 = np.zeros(len(vec))

            # first n values and last n values
            # first n values
            n = round(len(vec) / 10)
            y_sub = np.copy(vec[0:n])
            y_qual = np.zeros(len(y_sub))
            y_qual = rm_outliers.rm_IQR(y_sub, y_qual, factor=3)
            if len(idx1) != 0:
                y_qual0[0:n] = y_qual0[0:n] + y_qual  # update vec_qual

            # Last n values
            y_sub = np.copy(vec[-n:])
            y_qual = np.zeros(len(y_sub))
            y_qual = rm_outliers.rm_IQR(y_sub, y_qual, factor=3)
            if len(idx1) != 0:
                y_qual0[-n:] = y_qual0[-n:] + y_qual  # update vec_qual

            # visualization
            # filename = self.outfile[0:-3] + "_filter1.png"
            # rm_outliers.plotanomaly(vec, vec_qual, filename, Save=True)

            # Second filer:
            # remove points precedent filter
            ind = np.where(y_qual0 >= 1)[0]
            id_point = np.delete(id_point, ind)

            vec = np.copy(self.Temp)[id_point]
            y_qual1 = np.zeros(len(vec))

            if len(vec) >= 5:
                val_thresh = (np.quantile(vec, 0.9) * 2.5) - max(vec)
                y_qual1 = rm_outliers.rm_IQR_moving(vec, y_qual1, windowsize=round(len(vec) / 4), factor=10)
            else:
                val_thresh = -1
                #idx1 = []

            if val_thresh > 0 & len(idx1) != 0:
                id_point = np.delete(id_point, y_qual1)

            # Update vector PAR_qual
            self.PAR_qual = np.zeros(len(self.PAR_qual)) + 1
            self.PAR_qual[id_point] = 0

    def Temp_maintenance_extract():
        filename = "../renku_ppmooring/ppmooring/scripts/PAR_maintenance.csv"
        maintenance = pd.read_csv(filename, sep=";")

        vec = np.copy(maintenance.debut)
        min_interval = min(abs(np.diff(vec)))  # smallest time interval

        minvec = min(vec)
        maxvec = max(vec)

        vecnew = np.arange(minvec, maxvec, min_interval)

        # Find index of maintenance points
        xy, xind, yind = np.intersect1d(vec, vecnew, return_indices=True)

        # y=np.ones(len(vecnew))
        # plt.plot(np.arange(len(vecnew)),y)
        # plt.plot(yind,y[yind],'.k')

        # Create vec_index
        vec_index = np.zeros(len(vecnew))
        vec_index[yind] = 1

        # index of the end of each interval
        end = []
        for i in range(0, (len(vec_index) - 1)):
            val = vec_index[i] - vec_index[i + 1]
            if val == 1:
                end = np.append(end, i)

        end = end.astype(int)
        end = vecnew[end]

        # index of the start of each interval
        debut = np.array([0])
        for i in range(0, (len(vec_index) - 1)):
            val = vec_index[i] - vec_index[i + 1]
            if val == -1:
                debut = np.append(debut, i)
        debut = debut.astype(int)
        debut = vecnew[debut]

        # Create dataframe
        maintenance = {"debut": debut, "end": end}
        maintenance = pd.DataFrame(maintenance)

        return maintenance

    # Quality flags in Level 0
    def quality_flags_level0(self, save=False):
        ## Perform quality assurance tests on data
        ## Code quality: 0 = nothing to report, 1 = more investigation

        # negative values = 0:
        ind = np.where(self.Temp < 0)[0]
        if len(ind) != 0:
            self.Temp[ind] = 10000

        vec = np.copy(self.Temp)
        vec_qual = np.copy(self.Temp_qual)


        # Remove maintenance data
        # maintenance = DO_maintenance_extract()

        indanom = []  # index anomaly
        for a in range(0, len(self.maintenance["debut"])):
            ind0 = np.where(self.time >= self.maintenance["debut"][a])[0]
            ind1 = np.where(self.time <= self.maintenance["end"][a])[0]
            ind = np.intersect1d(ind0, ind1)

            if len(ind) != 0:
                indanom = np.append(indanom, ind)
                indanom = np.unique(indanom, axis=0)

        # First filter:
        # remove maintenance dates:
        if len(indanom) != 0:
            vec_qual[indanom.astype(int)] = 1

        # nan values
        ind = np.where(np.isnan(vec) == True)[0]
        if len(ind) != 0:
            vec[ind] = 10000

        # Remove values > 40 degC
        ind = np.where(vec >= 20)[0]
        if len(ind) != 0:
            vec_qual[ind] = 1

        # first n values and last n values
        # first n values
        n = 30
        y_sub = np.copy(vec[0:n])
        y_qual = np.zeros(len(y_sub))
        y_qual = rm_outliers.rm_IQR(y_sub, y_qual, factor=3)
        if sum(y_qual)>0:
            vec_qual[0:n][y_qual] = 1  # update vec_qual
            y_qual = []

        # Last n values
        y_sub = np.copy(vec[-n:])
        y_qual = np.zeros(len(y_sub))
        y_qual = rm_outliers.rm_IQR(y_sub, y_qual, factor=3)
        if sum(y_qual)>0:
            vec_qual[-n:][y_qual] = 1  # update vec_qual
            y_qual=[]

        # Remove threshold values error
        vdepth = self.depth[0]

        valid_range = {
            "50": [4, 40],
            "250": [4, 35],
            "500": [4, 35],
            "750": [4, 35],
            "1000": [4, 30],
            "1250": [4, 30],
            "1500": [4, 30],
            "1750": [4, 30],
            "2000": [4, 25],
            "2250": [4, 25],
            "2500": [4, 25],
            "2750": [4, 25],
            "3000": [4, 25],
            "5000": [2, 20],
            "10000": [2, 20],
        }

        # Check if variable is within valid range

        indlow = np.where(vec < valid_range[str(vdepth)][0])
        indhigh = np.where(vec >= valid_range[str(vdepth)][1])
        vec_qual[indlow] = 1
        vec_qual[indhigh] = 1

        # Update vector Temp_qual
        self.Temp_qual = np.copy(vec_qual)


    # Quality flags in Level 1A
    def quality_flags_level1A(self, save=False):
        ## Perform quality assurance tests on data
        ## Code quality: 0 = nothing to report, 1 = more investigation

        # Variable Temp
        vec = np.copy(self.Temp)
        vec_qual = np.copy(self.Temp_qual)
        id_point = np.array(range(0, len(vec)))

        # remove points precedent filter
        ind = np.where(vec_qual >= 1)[0]
        id_point = np.delete(id_point, ind)

        # Reinitialize
        vec = np.copy(self.Temp)[id_point]
        vec_qual = np.zeros(len(vec))

        if len(vec) < 10:
            # Final: all values superior than 1 in vec_qual become 1
            self.Temp_qual = np.zeros(len(self.Temp_qual)) + 1
            self.Temp_qual[id_point] = 0
            return

        # Next _filter
        xval = max(vec) / np.quantile(vec, 0.95)
        if xval > 1.5:
            y_qual = rm_outliers.rm_kmeans_threshold(vec, ncluster=2, threshold=5)
            if sum(y_qual)> 0:
                id_point = np.delete(id_point, y_qual)
                #idx = []
                # Reinitialize
                vec = np.copy(self.Temp)[id_point]
                vec_qual = np.zeros(len(vec))
                y_qual = []

        # Next _filter
        xval = max(vec) / np.quantile(vec, 0.95)
        if xval > 1.2:
            dayperiode = (24 * 60 * 60 / (self.time[1] - self.time[0])).astype(int)
            if len(vec) >= 3 * dayperiode:
                y_qual = rm_outliers.rm_convolution(vec, vec_qual, window_len=dayperiode, window_type='blackman',
                                                         threshold=np.quantile(vec, 0.9999), n_sigma=2.5)
            else:
                y_qual = rm_outliers.rm_convolution(vec, vec_qual, window_len=round(len(vec) / 2),
                                                         window_type='blackman',
                                                         threshold=np.quantile(vec, 0.9999), n_sigma=2.5)
            if sum(y_qual) > 0:
                id_point = np.delete(id_point, y_qual)
                #idx = []
                # Reinitialize
                vec = np.copy(self.Temp)[id_point]
                vec_qual = np.zeros(len(vec))
                y_qual=[]

        # Next _filter
        xval = max(vec) - np.quantile(vec, 0.95)
        if xval > 2:
            y_qual = np.where(vec >= np.quantile(vec, 0.9999))[0]

            if sum(y_qual) != 0:
                id_point = np.delete(id_point, y_qual)
                idx = []
                # Reinitialize
                vec = np.copy(self.Temp)[id_point]
                vec_qual = np.zeros(len(vec))
                y_qual = []

        xval = max(vec) - np.quantile(vec, 0.95)
        if xval > 2:
            y_qual = rm_outliers.rm_variation_rate(vec, vec_qual)

            if sum(y_qual) > 0:
                id_point = np.delete(id_point, y_qual)
                #idx = []
                # Reinitialize
                vec = np.copy(self.Temp)[id_point]
                vec_qual = np.zeros(len(vec))
                y_qual = []

        # Final: all values superior than 1 in vec_qual become 1
        self.Temp_qual = np.zeros(len(self.Temp_qual)) + 1
        self.Temp_qual[id_point] = 0

class rm_outliers(ppmooring):
    ## Class for removal outliers functions and visualization

    #visualization
    def plotanomaly(vec,vec_qual,filename="test.png",Save=False):
        plt.ioff()

        plt.scatter(np.arange(len(vec)), vec, c='b',label='raw')

        idx = np.where(vec_qual >= 1)[0]
        plt.scatter(idx, vec[idx], c='r', label='outlier')
        plt.title('Detected outliers')
        plt.xlabel('Time steps')
        plt.ylabel('Values')
        plt.legend();

        if Save:
            plt.savefig(filename)
        plt.close()

    def rm_IQR(vecX,vecX_qual,factor=3):
        """
            Remove values outside Interquartile Range (IQR) for timeseries data.

            Parameters:
                vecX (np.array): Data array to which to apply the quality assurance
                vecX_qual (np.array): An array of bools where True means non-trusted data before this filter
                factor (np.int): threshold for outlier labelling rule

            Returns:
                y_qual (np.array): An array of bools where True means non-trusted data after this outlier detection
        """
        y = np.copy(vecX)
        y_qual = np.copy(vecX_qual)
        q75 = np.quantile(y, 0.75)
        q25 = np.quantile(y, 0.25)
        IQR = q75 - q25
        if IQR != 0:
            outsup = q75 + factor * IQR
            outinf = q25 - factor * IQR

            idx1 = np.where(y >= outsup)[0]
            idx2 = np.where(y <= outinf)[0]

            idx0 = np.r_[idx1, idx2]
        else:
            idx0 = np.array([])

        if len(idx0) == 0:
            # print("No outliers detected - IQR")
            idx = idx0
        else:
            y_qual[idx0] = y_qual[idx0] + 1
            y_qual = np.array(y_qual,dtype=bool)

            #idx = np.where(y_qual >= 1)[0]
        return y_qual

    def rm_variation_rate(vec, vecX_qual):
        """
            Remove values if variation rate exceed a defined threshold.

            Parameters:
                vec (np.array): Data array to which to apply the quality assurance
                vecX_qual (np.array): An array of bools where True means non-trusted data

            Returns:
                y_qual (np.array): An array of bools where True means non-trusted data for this outlier dectection
        """
        y_qual = np.copy(vecX_qual)

        vecdiff = abs(np.diff(vec))
        if len(vecdiff)>0:
            vecdiff_quan = np.quantile(vecdiff,0.999)
            idx_vecdiff = np.where(vecdiff>=vecdiff_quan)[0]

            vec_max = np.max(vec)
            if vec_max / np.quantile(vec, 0.99) < 2:
                quantile_threshold = 0.99
            elif vec_max / np.quantile(vec, 0.999) < 2:
                quantile_threshold = 0.999
            elif vec_max / np.quantile(vec, 0.9999) < 2:
                quantile_threshold = 0.9999
            elif vec_max / np.quantile(vec, 0.99999) < 2:
                quantile_threshold = 0.99999

            vec_quan = np.quantile(vec, quantile_threshold)
            idx_vec = np.where(vec >= vec_quan)[0]
            idx = list(set(idx_vecdiff) & set(idx_vec))
        else:
            idx = []

        y_qual[idx] = y_qual[idx] + 1
        y_qual = np.array(y_qual, dtype=bool)

        return y_qual

    def rm_IQR_moving(vecX, vecX_qual, windowsize=15, factor=3, plot=False):
        """
           Remove outliers values based on Interquartile Range (IQR) for a window of time series data

           Parameters:
               vec (np.array): Data array to which to apply the quality assurance
               vecX_qual (np.array): An array of bools where True means non-trusted data
               windowsize (np.int): window size of data

           Returns:
               y_qual (np.array): An array of bools where True means non-trusted data for this outlier dectection
       """
        y = np.copy(vecX)
        y_qual = np.copy(vecX_qual)
        N_y = len(y)

        if N_y < windowsize:
            print("window size is larger than len(y)")
        else:
            for i in range(0, (N_y - windowsize + 1)):
                y_sub = np.copy(y[i:i + windowsize])

                q75 = np.quantile(y_sub, 0.75)
                q25 = np.quantile(y_sub, 0.25)
                IQR = q75 - q25
                outsup = q75 + factor * IQR
                outinf = q25 - factor * IQR

                idx1 = np.where(y_sub >= outsup)[0]
                idx2 = np.where(y_sub <= outinf)[0]

                idx0 = np.r_[idx1, idx2]

                if len(idx0) != 0:
                    y_qual[i + idx0] = y_qual[i + idx0] + 1

        y_qual = np.array(y_qual, dtype=bool)

        return y_qual

    def rm_max(vecX, vecX_qual, factor1=3, semiwindow=1000):
        """
           Remove outliers values based on Interquartile Range (IQR) for a window of time series data

           Parameters:
               vec (np.array): Data array to which to apply the quality assurance
               vecX_qual (np.array): An array of bools where True means non-trusted data
               semiwindow (int): window size of data
               factor1 (int): threshold for outlier labelling rule

           Returns:
               y_qual (np.array): An array of bools where True means non-trusted data for this outlier dectection
       """

        y = np.copy(vecX)
        y_qual = np.copy(vecX_qual)

        maxy = np.nanmax(y)

        n0 = np.where(y == maxy)[0][0] - semiwindow
        n1 = np.where(y == maxy)[0][0] + semiwindow

        if n0 < 0:
            n0 = 0
        if n1 > (len(y) - 1):
            n1 = len(y) - 1

        vec_sub = y[n0:n1]
        vec_qual = np.zeros(len(vec_sub))

        y_qual99 = rm_outliers.rm_IQR(vec_sub, vec_qual, factor=factor1)
        if sum(y_qual99)>0:
            y_qual[n0:n1] = y_qual[n0:n1] + y_qual99  # update vec_qual

        y_qual = np.array(y_qual, dtype=bool)
        return y_qual

    def rm_convolution(vecX,vecX_qual,window_len = 30,window_type = 'blackman',n_sigma =2,threshold=20):
        """
            Remove values using convolutional smoothing of single or multiple time-series

            Parameters:
                vecX (np.array): Data array to which to apply the quality assurance
                vecX_qual (np.array): An array of bools where True means non-trusted data
                window_len (int) : Greater than equal to 1. The length of the window used to compute
        the convolutions.
                window_type (str):  The type of the window used to compute the convolutions.
        Supported types are: 'ones', 'hanning', 'hamming', 'bartlett', 'blackman'.

            Returns:
                y_qual (np.array): An array of bools where True means non-trusted data for this outlier detection
        """

        timesteps = len(vecX)

        series = defaultdict(partial(np.ndarray, shape=(1), dtype='float32'))

        for i in tqdm(range(timesteps + 1), total=(timesteps + 1)):
            if i > window_len:
                smoother = ConvolutionSmoother(window_len=window_len, window_type=window_type)
                smoother.smooth(series['original'][-window_len:])

                series['smooth'] = np.hstack([series['smooth'], smoother.smooth_data[:, -1]])

                _low, _up = smoother.get_intervals('sigma_interval', n_sigma=n_sigma)
                series['low'] = np.hstack([series['low'], _low[:, -1]])
                series['up'] = np.hstack([series['up'], _up[:, -1]])

                is_anomaly = np.logical_or(
                    series['original'][-1] > series['up'][-1],
                    series['original'][-1] < series['low'][-1]
                )

                if is_anomaly.any():
                    series['idx'] = np.hstack([series['idx'], is_anomaly * i]).astype(int)

            if i >= timesteps:
                continue

            series['original'] = np.hstack([series['original'], vecX[i]])

        if len(series["idx"])!=0:
            idx0 = np.where(series['original']>threshold)[0]
            idx =  np.intersect1d(idx0,series['idx'])
            if len(idx)!=0:
                if idx[-1] == len(vecX):
                    idx[-1] = idx[-1]-1
        y_qual = np.copy(vecX_qual)

        if len(idx) !=0:
            y_qual[idx] = 1
            y_qual = np.array(y_qual,dtype=bool)

        return y_qual

    def rm_kmeans(vecX, ncluster=2):
        """
            Remove outliers based on kmean clustering.

            Parameters:
                vecX (np.array): Data array to which to apply the quality assurance
                ncluster (int) : number of cluster (>=2)
            Returns:
                y_qual (np.array): An array of bools where True means non-trusted data for this outlier detection
        """
        y_qual = np.zeros(len(vecX), dtype=bool)
        clusterer = KMeans(n_clusters=ncluster)
        clusterer.fit(vecX.reshape(-1, 1))

        nearest_centroid_idx = clusterer.predict(vecX.reshape(-1, 1))

        igr1 = np.where(nearest_centroid_idx == 0)[0]
        igr2 = np.where(nearest_centroid_idx == 1)[0]

        val_thresh = (np.mean(vecX[igr2]) - np.mean(vecX[igr1])) / np.quantile(vecX, 0.90)

        if val_thresh >= 5:  # if there is no 2 clearly seperated groups
            y_qual[igr2] = True

        y_qual = np.array(y_qual, dtype=bool)

        return y_qual

    def rm_kmeans_threshold(vecX, ncluster=2,threshold = 1.2):
        """
            Remove outliers based on kmean clustering and threshold value.

            Parameters:
                vecX (np.array): Data array to which to apply the quality assurance
                ncluster : number of cluster (>=2)
            Returns:
                y_qual (np.array): An array of bools where True means non-trusted data for this outlier detection
        """
        y_qual = np.zeros(len(vecX), dtype=bool)
        clusterer = KMeans(n_clusters=ncluster)
        clusterer.fit(vecX.reshape(-1, 1))

        nearest_centroid_idx = clusterer.predict(vecX.reshape(-1, 1))

        igr1 = np.where(nearest_centroid_idx == 0)[0]
        igr2 = np.where(nearest_centroid_idx == 1)[0]

        if len(igr1)>len(igr2):
            igrfin = igr2
        else:
            igrfin = igr1

        val_thresh = abs((np.mean(vecX[igr2]) - np.mean(vecX[igr1]))) / np.quantile(vecX, 0.90)

        if val_thresh < threshold:  # if there is no 2 clearly seperated groups
            igrfin = []
        else:
            y_qual[igrfin] = True

        return y_qual

