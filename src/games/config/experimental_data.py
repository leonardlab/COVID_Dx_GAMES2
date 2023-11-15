#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 10:36:16 2022

@author: kate
"""
from typing import Tuple, List
import pandas as pd
import pickle


def define_experimental_data(settings: dict) -> Tuple[List[float], List[float], List[float]]:
   """
   Imports experimental data

   Parameters
   ---------
   settings
      a dictionary of run settings

   Returns
   -------
   x
      a list of floats defining the independent variable for the given dataset

   exp_data
      a list of floats defining the normalized dependent variable for the given
      dataset

   exp_error
      a list of floats defining the normalized measurement error for
      the dependent variable for the given dataset
   """

   path = settings["context"] + "config/"
   # filename = path + "training_data_" + settings["dataID"] + ".csv"
   # df_exp = pd.read_csv(filename)
   # x = list(df_exp["x"])
   # exp_data = list(df_exp["y"])
   # exp_error = list(df_exp["y_err"])

   if "rep1" in settings["dataID"]:
      filename_data = path + "PROCESSED DATA_EXP.pkl"
      filename_error = path + "PROCESSED DATA_ERR.pkl"

   elif "rep2" in settings["dataID"]:
      filename_data = path + "PROCESSED_DATA_rep2_EXP.pkl"
      filename_error = path + "PROCESSED_DATA_rep2_ERR.pkl"

   if "rep3" in settings["dataID"]:
      filename_data = path + "PROCESSED_DATA_rep3_EXP.pkl"
      filename_error = path + "PROCESSED_DATA_rep3_ERR.pkl"

   df_data = pd.read_pickle(filename_data)
   df_error = pd.read_pickle(filename_error)

#Choose conditions to include or drop                             
   if settings["dataID"] == ' rep1 all echo drop high error':
      drop_labels = [
         [5.0, 10.0, 0.001, 1, 90], [1.0, 2.5, 0.001, 10.0, 90.0], [20.0, 10.0, 0.001, 1.0, 90.0],
         [5.0, 0.5, 0.005, 10.0, 90.0], [20.0, 0.5, 0.005, 1.0, 90.0], [1.0, 2.5, 0.001, 1.0, 90.0], 
         [20.0, 0.5, 0.005, 10.0, 90.0]
      ]
      
      x = []
      exp_data = []
      error = []
      timecourses = []
      timecourses_err = []
      maxVal = 0.6599948235700113
      
      for (columnName, columnData) in df_data.iteritems():
         label = list(columnData.iloc[0])
         if (label == [20.0, 10.0, 0.001, 10.0, 90.0] or label == [20.0, 10.0, 0.001, 0.0, 90.0]
             or label == [5.0, 2.5, 0.001, 0.0, 90.0]):
               continue
         
         elif label in drop_labels:
               continue
      
         else:
            
               x.append(label)
         
               timecourse = list(columnData.iloc[1:])
               timecourse = [i/maxVal for i in timecourse]
               timecourses.append(timecourse)
               exp_data = exp_data + timecourse
      
      for (columnName, columnData) in df_error.iteritems():
         label = list(columnData.iloc[0])
         if label == ([20.0, 10.0, 0.001, 10.0, 90.0] or label == [20.0, 10.0, 0.001, 0.0, 90.0]
                      or label == [5.0, 2.5, 0.001, 0.0, 90.0]):
               continue
         
         elif label in drop_labels:
               continue
         
         else:
               err = list(columnData.iloc[1:])
               err = [i/maxVal for i in err]
               timecourses_err.append(err)
               error = error + err 

   elif settings["dataID"] == 'rep1 all echo not in slice drop high error':
      labels1 = [[1.0, 2.5, 0.005, 1, 90], [5.0, 2.5, 0.005, 1, 90], [20.0, 2.5, 0.005, 1, 90]]
      labels10 = [[1.0, 2.5, 0.005, 10, 90], [5.0, 2.5, 0.005, 10, 90], [20.0, 2.5, 0.005, 10, 90]]
      labels_T7 = labels1 + labels10
      
      labels1 = [[5.0, 0.5, 0.005, 1, 90], [5.0, 2.5, 0.005, 1, 90], [5.0, 10.0, 0.005, 1, 90]]
      labels10 = [[5.0, 0.5, 0.005, 10, 90], [5.0, 2.5, 0.005, 10, 90], [5.0, 10.0, 0.005, 10, 90]]
      labels_RT = labels1 + labels10
      
      labels1 = [[5.0, 2.5, 0.001, 1, 90], [5.0, 2.5, 0.005, 1, 90], [5.0, 2.5, 0.02, 1, 90]]
      labels10 = [[5.0, 2.5, 0.001, 10, 90], [5.0, 2.5, 0.005, 10, 90], [5.0, 2.5, 0.02, 10, 90]]
      labels_RNase =labels1 + labels10
      
      labels_slice = labels_T7 + labels_RT + labels_RNase

      drop_labels = [
         [5.0, 10.0, 0.001, 1, 90], [1.0, 2.5, 0.001, 10.0, 90.0], [20.0, 10.0, 0.001, 1.0, 90.0],
         [5.0, 0.5, 0.005, 10.0, 90.0], [20.0, 0.5, 0.005, 1.0, 90.0], [1.0, 2.5, 0.001, 1.0, 90.0], 
         [20.0, 0.5, 0.005, 10.0, 90.0]
      ]

      x = []
      exp_data = []
      error = []
      timecourses = []
      timecourses_err = []
      maxVal = 0.6599948235700113

      for (columnName, columnData) in df_data.iteritems():
         label = list(columnData.iloc[0])
         if label == ([20.0, 10.0, 0.001, 10.0, 90.0] or label == [20.0, 10.0, 0.001, 0.0, 90.0]
                      or label == [5.0, 2.5, 0.001, 0.0, 90.0]):
               continue

         elif label in labels_slice or label in drop_labels:
               continue

         else:
               x.append(label)
               timecourse = list(columnData.iloc[1:])
               timecourse = [i/maxVal for i in timecourse]
               timecourses.append(timecourse)
               exp_data = exp_data + timecourse

      for (columnName, columnData) in df_error.iteritems():
         label = list(columnData.iloc[0])
         if label == ([20.0, 10.0, 0.001, 10.0, 90.0] or label == [20.0, 10.0, 0.001, 0.0, 90.0]
                      or label == [5.0, 2.5, 0.001, 0.0, 90.0]):
               continue

         elif label in labels_slice or label in drop_labels:
               continue

         else:
               err = list(columnData.iloc[1:])
               err = [i/maxVal for i in err]
               timecourses_err.append(err)
               error = error + err 


   elif settings["dataID"] == 'rep 1 all echo':
      x = []
      exp_data = []
      error = []
      timecourses = []
      timecourses_err = []
      maxVal = 0.6599948235700113
      
      for (columnName, columnData) in df_data.iteritems():
         label = list(columnData.iloc[0])
         if label == ([20.0, 10.0, 0.001, 10.0, 90.0] or label == [20.0, 10.0, 0.001, 0.0, 90.0]
                      or label == [5.0, 2.5, 0.001, 0.0, 90.0]):
               continue
      
         else:
               x.append(label)
               timecourse = list(columnData.iloc[1:])
               timecourse = [i/maxVal for i in timecourse]
               timecourses.append(timecourse)
               exp_data = exp_data + timecourse
      
      for (columnName, columnData) in df_error.iteritems():
         label = list(columnData.iloc[0])
         if label == ([20.0, 10.0, 0.001, 10.0, 90.0] or label == [20.0, 10.0, 0.001, 0.0, 90.0]
                      or label == [5.0, 2.5, 0.001, 0.0, 90.0]):
               continue
         else:
               err = list(columnData.iloc[1:])
               err = [i/maxVal for i in err]
               timecourses_err.append(err)
               error = error + err   
                
    
   elif settings["dataID"] == 'rep 1 slice drop high error':
      labels1 = [[1.0, 2.5, 0.005, 1, 90], [5.0, 2.5, 0.005, 1, 90], [20.0, 2.5, 0.005, 1, 90]]
      labels10 = [[1.0, 2.5, 0.005, 10, 90], [5.0, 2.5, 0.005, 10, 90], [20.0, 2.5, 0.005, 10, 90]]
      labels_T7 = labels1 + labels10
      
      labels1 = [[5.0, 0.5, 0.005, 1, 90], [5.0, 2.5, 0.005, 1, 90], [5.0, 10.0, 0.005, 1, 90]]
      labels10 = [[5.0, 0.5, 0.005, 10, 90], [5.0, 2.5, 0.005, 10, 90], [5.0, 10.0, 0.005, 10, 90]]
      labels_RT = labels1 + labels10
      
      labels1 = [[5.0, 2.5, 0.001, 1, 90], [5.0, 2.5, 0.005, 1, 90], [5.0, 2.5, 0.02, 1, 90]]
      labels10 = [[5.0, 2.5, 0.001, 10, 90], [5.0, 2.5, 0.005, 10, 90], [5.0, 2.5, 0.02, 10, 90]]
      labels_RNase =labels1 + labels10
      
      labels_pre_drop = labels_T7 + labels_RT + labels_RNase
      drop_labels = [
         [5.0, 10.0, 0.001, 1, 90], [1.0, 2.5, 0.001, 10.0, 90.0], [20.0, 10.0, 0.001, 1.0, 90.0],
         [5.0, 0.5, 0.005, 10.0, 90.0], [20.0, 0.5, 0.005, 1.0, 90.0], [1.0, 2.5, 0.001, 1.0, 90.0],
         [20.0, 0.5, 0.005, 10.0, 90.0]
      ]
      labels = []
      for label in labels_pre_drop:
         if label not in drop_labels:
               labels.append(label)
   
      x = []
      exp_data = []
      for (columnName, columnData) in df_data.iteritems():
         label = list(columnData.iloc[0])
         
         if label == ([20.0, 10.0, 0.001, 10.0, 90.0] or label == [20.0, 10.0, 0.001, 0.0, 90.0]
                      or label == [5.0, 2.5, 0.001, 0.0, 90.0]):
               continue
         elif label in labels:
               x.append(label)
               timecourse = list(columnData.iloc[1:])

               exp_data = exp_data + timecourse
      
      maxVal = max(exp_data)
      exp_data = []
      timecourses = []
      timecourses_err = []
      for (columnName, columnData) in df_data.iteritems():
         label = list(columnData.iloc[0])
         
         if label == ([20.0, 10.0, 0.001, 10.0, 90.0] or label == [20.0, 10.0, 0.001, 0.0, 90.0]
                      or label == [5.0, 2.5, 0.001, 0.0, 90.0]):
               continue
         elif label in labels:
               timecourse = list(columnData.iloc[1:])
               timecourse = [i/maxVal for i in timecourse]
               
               timecourses.append(timecourse)
               exp_data = exp_data + timecourse           
   
      error = []
      for (columnName, columnData) in df_error.iteritems():
         label = list(columnData.iloc[0])
         if label == ([20.0, 10.0, 0.001, 10.0, 90.0] or label == [20.0, 10.0, 0.001, 0.0, 90.0]
                      or label == [5.0, 2.5, 0.001, 0.0, 90.0]):
               continue
         
         elif label in labels:
               err = list(columnData.iloc[1:])
               err = [i/maxVal for i in err]
               timecourses_err.append(err)
               error = error + err 


   elif settings["dataID"] == 'rep2 all echo drop high error':
      drop_labels = [[20.0, 2.5, 0.02, 10, 90]]        
      x = []
      exp_data = []
      error = []
      timecourses = []
      timecourses_err = []
      maxVal = 2.94995531724754
      
      for (columnName, columnData) in df_data.iteritems():
         label = list(columnData.iloc[0])
         
         if label in drop_labels:
               continue
         else:
            
               x.append(label)
         
               timecourse = list(columnData.iloc[1:])
               timecourse = [i/maxVal for i in timecourse]
               timecourses.append(timecourse)
               exp_data = exp_data + timecourse
      
      for (columnName, columnData) in df_error.iteritems():
         label = list(columnData.iloc[0])
         
         if label in drop_labels:
               continue
         
         else:
               err = list(columnData.iloc[1:])
               err = [i/maxVal for i in err]
               timecourses_err.append(err)
               error = error + err    


   elif settings["dataID"] == 'rep2 all echo not in slice drop high error':
      labels1 = [[1.0, 2.5, 0.005, 1, 90], [5.0, 2.5, 0.005, 1, 90], [20.0, 2.5, 0.005, 1, 90]]
      labels10 = [[1.0, 2.5, 0.005, 10, 90], [5.0, 2.5, 0.005, 10, 90], [20.0, 2.5, 0.005, 10, 90]]
      labels_T7 = labels1 + labels10
      
      labels1 = [[5.0, 0.5, 0.005, 1, 90], [5.0, 2.5, 0.005, 1, 90], [5.0, 10.0, 0.005, 1, 90]]
      labels10 = [[5.0, 0.5, 0.005, 10, 90], [5.0, 2.5, 0.005, 10, 90], [5.0, 10.0, 0.005, 10, 90]]
      labels_RT = labels1 + labels10
      
      labels1 = [[5.0, 2.5, 0.001, 1, 90], [5.0, 2.5, 0.005, 1, 90], [5.0, 2.5, 0.02, 1, 90]]
      labels10 = [[5.0, 2.5, 0.001, 10, 90], [5.0, 2.5, 0.005, 10, 90], [5.0, 2.5, 0.02, 10, 90]]
      labels_RNase =labels1 + labels10
      
      labels_slice = labels_T7 + labels_RT + labels_RNase

      drop_labels = [[20.0, 2.5, 0.02, 10, 90]]

      x = []
      exp_data = []
      error = []
      timecourses = []
      timecourses_err = []
      maxVal = 2.94995531724754

      for (columnName, columnData) in df_data.iteritems():
         label = list(columnData.iloc[0]) 
         if label in labels_slice or label in drop_labels:
               continue

         else:
               x.append(label)
               timecourse = list(columnData.iloc[1:])
               timecourse = [i/maxVal for i in timecourse]
               timecourses.append(timecourse)
               exp_data = exp_data + timecourse

      for (columnName, columnData) in df_error.iteritems():
         label = list(columnData.iloc[0])
         if label in labels_slice or label in drop_labels:
               continue

         else:
               err = list(columnData.iloc[1:])
               err = [i/maxVal for i in err]
               timecourses_err.append(err)
               error = error + err 


   elif settings["dataID"] == 'rep2 slice drop high error':
      labels1 = [[1.0, 2.5, 0.005, 1, 90], [5.0, 2.5, 0.005, 1, 90], [20.0, 2.5, 0.005, 1, 90]]
      labels10 = [[1.0, 2.5, 0.005, 10, 90], [5.0, 2.5, 0.005, 10, 90], [20.0, 2.5, 0.005, 10, 90]]
      labels_T7 = labels1 + labels10
      
      labels1 = [[5.0, 0.5, 0.005, 1, 90], [5.0, 2.5, 0.005, 1, 90], [5.0, 10.0, 0.005, 1, 90]]
      labels10 = [[5.0, 0.5, 0.005, 10, 90], [5.0, 2.5, 0.005, 10, 90], [5.0, 10.0, 0.005, 10, 90]]
      labels_RT = labels1 + labels10
      
      labels1 = [[5.0, 2.5, 0.001, 1, 90], [5.0, 2.5, 0.005, 1, 90], [5.0, 2.5, 0.02, 1, 90]]
      labels10 = [[5.0, 2.5, 0.001, 10, 90], [5.0, 2.5, 0.005, 10, 90], [5.0, 2.5, 0.02, 10, 90]]
      labels_RNase =labels1 + labels10
      
      labels_pre_drop = labels_T7 + labels_RT + labels_RNase
      drop_labels = [[20.0, 2.5, 0.02, 10, 90]] #not in the slice labels so not relevant here
      labels = []
      for label in labels_pre_drop:
         if label not in drop_labels:
               labels.append(label)
   
      x = []
      exp_data = []
      for (columnName, columnData) in df_data.iteritems():
         label = list(columnData.iloc[0])
      
         if label in labels:
               x.append(label)
               timecourse = list(columnData.iloc[1:])

               exp_data = exp_data + timecourse
      
      maxVal = max(exp_data)
      exp_data = []
      timecourses = []
      timecourses_err = []
      for (columnName, columnData) in df_data.iteritems():
         label = list(columnData.iloc[0])

         if label in labels:
               timecourse = list(columnData.iloc[1:])
               timecourse = [i/maxVal for i in timecourse]
               
               timecourses.append(timecourse)
               exp_data = exp_data + timecourse           
   
      error = []
      for (columnName, columnData) in df_error.iteritems():
         label = list(columnData.iloc[0])
         
         if label in labels:
               err = list(columnData.iloc[1:])
               err = [i/maxVal for i in err]
               timecourses_err.append(err)
               error = error + err 

   elif settings["dataID"] == 'rep3 all echo drop high error':
      drop_labels = [[5.0, 10.0, 0.02, 10, 90]]        
      x = []
      exp_data = []
      error = []
      timecourses = []
      timecourses_err = []
      maxVal = 1.12314566577301
      
      for (columnName, columnData) in df_data.iteritems():
         label = list(columnData.iloc[0])
         
         if label in drop_labels:
               continue
         else:
            
               x.append(label)
         
               timecourse = list(columnData.iloc[1:])
               timecourse = [i/maxVal for i in timecourse]
               timecourses.append(timecourse)
               exp_data = exp_data + timecourse
      
      for (columnName, columnData) in df_error.iteritems():
         label = list(columnData.iloc[0])
         
         if label in drop_labels:
               continue
         
         else:
               err = list(columnData.iloc[1:])
               err = [i/maxVal for i in err]
               timecourses_err.append(err)
               error = error + err    


   elif settings["dataID"] == 'rep3 all echo not in slice drop high error':
      labels1 = [[1.0, 2.5, 0.005, 1, 90], [5.0, 2.5, 0.005, 1, 90], [20.0, 2.5, 0.005, 1, 90]]
      labels10 = [[1.0, 2.5, 0.005, 10, 90], [5.0, 2.5, 0.005, 10, 90], [20.0, 2.5, 0.005, 10, 90]]
      labels_T7 = labels1 + labels10
      
      labels1 = [[5.0, 0.5, 0.005, 1, 90], [5.0, 2.5, 0.005, 1, 90], [5.0, 10.0, 0.005, 1, 90]]
      labels10 = [[5.0, 0.5, 0.005, 10, 90], [5.0, 2.5, 0.005, 10, 90], [5.0, 10.0, 0.005, 10, 90]]
      labels_RT = labels1 + labels10
      
      labels1 = [[5.0, 2.5, 0.001, 1, 90], [5.0, 2.5, 0.005, 1, 90], [5.0, 2.5, 0.02, 1, 90]]
      labels10 = [[5.0, 2.5, 0.001, 10, 90], [5.0, 2.5, 0.005, 10, 90], [5.0, 2.5, 0.02, 10, 90]]
      labels_RNase =labels1 + labels10
      
      labels_slice = labels_T7 + labels_RT + labels_RNase

      drop_labels = [[5.0, 10.0, 0.02, 10, 90]]

      x = []
      exp_data = []
      error = []
      timecourses = []
      timecourses_err = []
      maxVal = 1.12314566577301

      for (columnName, columnData) in df_data.iteritems():
         label = list(columnData.iloc[0]) 
         if label in labels_slice or label in drop_labels:
               continue

         else:
               x.append(label)
               timecourse = list(columnData.iloc[1:])
               timecourse = [i/maxVal for i in timecourse]
               timecourses.append(timecourse)
               exp_data = exp_data + timecourse

      for (columnName, columnData) in df_error.iteritems():
         label = list(columnData.iloc[0])
         if label in labels_slice or label in drop_labels:
               continue

         else:
               err = list(columnData.iloc[1:])
               err = [i/maxVal for i in err]
               timecourses_err.append(err)
               error = error + err 


   elif settings["dataID"] == 'rep3 slice drop high error':
      labels1 = [[1.0, 2.5, 0.005, 1, 90], [5.0, 2.5, 0.005, 1, 90], [20.0, 2.5, 0.005, 1, 90]]
      labels10 = [[1.0, 2.5, 0.005, 10, 90], [5.0, 2.5, 0.005, 10, 90], [20.0, 2.5, 0.005, 10, 90]]
      labels_T7 = labels1 + labels10
      
      labels1 = [[5.0, 0.5, 0.005, 1, 90], [5.0, 2.5, 0.005, 1, 90], [5.0, 10.0, 0.005, 1, 90]]
      labels10 = [[5.0, 0.5, 0.005, 10, 90], [5.0, 2.5, 0.005, 10, 90], [5.0, 10.0, 0.005, 10, 90]]
      labels_RT = labels1 + labels10
      
      labels1 = [[5.0, 2.5, 0.001, 1, 90], [5.0, 2.5, 0.005, 1, 90], [5.0, 2.5, 0.02, 1, 90]]
      labels10 = [[5.0, 2.5, 0.001, 10, 90], [5.0, 2.5, 0.005, 10, 90], [5.0, 2.5, 0.02, 10, 90]]
      labels_RNase =labels1 + labels10
      
      labels_pre_drop = labels_T7 + labels_RT + labels_RNase
      drop_labels = [[5.0, 10.0, 0.02, 10, 90]] #not in the slice labels so not relevant here
      labels = []
      for label in labels_pre_drop:
         if label not in drop_labels:
               labels.append(label)
     
      x = []
      exp_data = []
      for (columnName, columnData) in df_data.iteritems():
         label = list(columnData.iloc[0])
      
         if label in labels:
               x.append(label)
               timecourse = list(columnData.iloc[1:])

               exp_data = exp_data + timecourse
      
      maxVal = max(exp_data)
      exp_data = []
      timecourses = []
      timecourses_err = []
      for (columnName, columnData) in df_data.iteritems():
         label = list(columnData.iloc[0])

         if label in labels:
               timecourse = list(columnData.iloc[1:])
               timecourse = [i/maxVal for i in timecourse]
               
               timecourses.append(timecourse)
               exp_data = exp_data + timecourse           
   
      error = []
      for (columnName, columnData) in df_error.iteritems():
         label = list(columnData.iloc[0])
         
         if label in labels:
               err = list(columnData.iloc[1:])
               err = [i/maxVal for i in err]
               timecourses_err.append(err)
               error = error + err

   return x, exp_data, error
