#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 15:17:53 2022

@author: kate
"""
from typing import List, Tuple, Generator
import numpy as np
import pandas as pd
from lmfit import Parameters, minimize
import matplotlib.pyplot as plt
from games.utilities.metrics import calc_r_sq

def plotModelingObjectives123(solutions: list) -> None:
    """Plots the modeling objectives involving Hill fit summary metrics
        (1, 2, 3). Can also be used to plot the Hill fit summary metrics 
        for the experimental data.
    
    Parameters
    ---------- 
    solutions
        a list of floats defining the data for each condition and 
        timepoint (length = # data points total, for all doses at 
        all timepoints)

    Returns
    -------
    None

    """
    def fitHill(y_exp: list) -> Tuple[float, float, float, float, float]:
        """Fits data to a hill function for a single set of
            conditions (component doses).

        Parameters
        ---------
        y_exp
            a list floats defining the normalized simulation 
            values for a single set of conditions

        Returns
        -------
        f0
            a float defining the initial value in the dataset

        fmax
            a float defining the final value in the dataset

        km
            a float defining the optimized parameter for t1/2

        n
            a float defining the optimized parameter for the 
            Hill coefficient

        R_sq
            a float defining the R squared value between the 
            data and the Hill fit

        """
        x = list(np.linspace(0, 240, 61)) #time (min)
        
        #Set v max to the final value of the time course
        fmax = y_exp[-1]
       
        #Set v0 to the intiial value of the time course
        f0 = y_exp[0]
    
        def residual(p: list, x: list, y_exp: list) -> float:
            """Calculates the residual between the input simulation
                values and the Hill fit. Used in the minimization
                function as the cost function to be minimized
                between the simulation and the Hill fit.

            Parameters
            ----------
            p
                a list of floats defining the parameters for
                the hill function

            x
                a list of floats defining the time values for
                the simulation

            y_exp
                a list of floats defining the simulation
                values

            Returns
            -------
            (y_exp - model)
                a float defining the residual between the
                input simulation values and the Hill fit

            """
            km = p['km'].value
            n = p['n'].value
            model = (((fmax - f0) * (x ** n)) / (km ** n + x ** n)) + f0
            return (y_exp - model)
        
        #Define parameters to be fit, along with initial guesses
        p = Parameters()
        p.add('km', value=10e-2, min=0, vary = True)
        p.add('n', value=2, min=0, max = 5, vary = True)
        
        #Perform the fit
        out = minimize(residual, p, args=(x, y_exp))
        bestFitParams = out.params.valuesdict()
        bestFitParamsList= bestFitParams.values()
        
        #Define fit parameters
        fit = []
        for value in bestFitParamsList:
            fit.append(value)
        km = fit[0]
        n = fit[1]
        
        #Simulate the Hill fit
        y_Hill = []
        for item in x:
            value = (((fmax - f0) * (item ** n)) / (km ** n + item ** n)) + f0
            y_Hill.append(value)
        
        #Calculate the R2 between the data and the Hill fit
        R_sq = calc_r_sq(y_exp, y_Hill)
        
        return f0, fmax, km, n, R_sq

    #Unnpack the data from "exp_data" and "error"
    def chunks(lst: list, n: int) -> Generator[list]:
        """Yield successive n-sized chunks from lst
        
        Parameters
        ----------
        lst
            a list of floats
            
        n
            an integer defining the size of each chunk

        Returns
        -------
        lst[i:i + n]
            a list of lists containing the values
            from lst structured as n-sized chunks

        """
        for i in range(0, len(lst), n):
            yield lst[i:i + n]
        
    def fitHillAllData(data_lists: list) -> None:
        """Fits data to a hill function for all conditions
            (component doses).

        Parameters
        ----------
        data_lists
            a list of lists containing the simulation values
            for each set of conditions
        
        Returns
        -------
        None

        Figures
        -------
        ''MODELING OBJECTIVES 123 ' + '.svg''
            Plot of modeling objectives involving Hill fit
            summary metrics (1, 2, 3)

        """
        fit_params = []
        for time_course in data_lists:
            f0, fmax, km, n, R_sq = fitHill(time_course)
            fit_params.append([f0, fmax, km, n, R_sq])
    
        df = pd.DataFrame()
        f0 = [list_[0] for list_ in fit_params]
        fmax = [list_[1] for list_ in fit_params]
        km = [list_[2] for list_ in fit_params]
        n = [list_[3] for list_ in fit_params]
        R_sq = [list_[4] for list_ in fit_params]
    
        df['f0'] = f0
        df['fmax'] = fmax
        df['km'] = km
        df['n'] = n
        df['R_sq'] = R_sq
        
        filename = 'experimental_summary_metrics'
        with pd.ExcelWriter(filename + '.xlsx') as writer:  
            df.to_excel(writer, sheet_name=' ')
    
        fig = plt.figure(figsize = (8,5))
        ax1 = plt.subplot(231)   
        ax2 = plt.subplot(232)
        ax3 = plt.subplot(233)
        ax4 = plt.subplot(234)
        ax5 = plt.subplot(235)
   
        color_ = 'dimgrey'
        bins_ = 25
        ax1.hist(f0, bins=bins_, color = color_)
        ax1.set_xlabel('f0')
    
        ax2.hist(fmax, bins=bins_, color = color_)
        ax2.set_xlabel('fmax')
      
        ax3.hist(km, bins=bins_, color = color_)
        ax3.set_xlabel('t 1/2')
        #ax3.set_xlim([50, 90])
    
        ax4.hist(n, bins=bins_, color = color_)
        ax4.set_xlabel('n')
    
        ax5.hist(R_sq, bins=bins_, color = color_)
        ax5.set_xlabel('R_sq')
        #ax5.set_xlim([0.994, 0.999])
    
        plt.savefig('modeling_objectives_123' + '.svg', dpi = 600)
 
    data_lists = list(chunks(solutions, 61))
    print('Determining Hill fits...')
    fitHillAllData(data_lists)
    print('Modeling objectives 1 2 and 3 plotted.')

def resultsPanel(
        dfSim: pd.DataFrame, dfExp: pd.DataFrame,
        dfErr: pd.DataFrame, dataID: str, labels: list,
        varyCondition: str, maxVal: float, y_max_RT: list, 
        y_max_T7: list, y_max_RNase: list
 ) -> None:   
    """Plots selected readout time courses (sweeps) for a given enzyme
    
    Parameters
    ---------- 
    dfSim
        a dataframe containing the simulated data

    dfExp
        a dataframe containing the experimental data

    dfErr
        a dataframe containing the measurement error associated with
        the experimental data
    
    labels
        a list of lists containing the condition labels to be plotted

    varyCondition
        a string contraining the label of the enzyme condition
        that is sweeped over in the plot
   
    Returns
    -------
    None
        
    Figures
    ------- 
    'modeling_objective_' + str(objective) + '.svg''
        plot of readout dynamics associated with the given modeling objective
        (enzyme sweep)

    """
    fig = plt.figure(figsize = (12,3))
    fig.subplots_adjust(wspace=0.1, hspace=0.1)
    ax1 = plt.subplot(141)   
    ax2 = plt.subplot(142)
    ax3 = plt.subplot(143)   
    ax4 = plt.subplot(144)
   
    def grabData(labels: list) -> Tuple[list, list, list]:
        """Compiles the simulation values and experimental
            data for the given set of labels

        Parameters
        ----------
        labels
            a list of lists defining the labels
            (component doses) to compile data for
            
        Returns
        -------
        sim
            a list of lists defining the simulation
            values for the set of labels

        exp
            a list of lists defining the experimental
            data for the set of labels

        err
            a list of lists defining the experimental
            error for the set of labels

        """
        count = 0
        sim = [[0] * 61] * 3
           
        for (columnName, columnData) in dfSim.items():
            label = str(columnName)
    
            time = np.linspace(0, 240, 61)
            t_course = columnData
   
            if label == str(labels[0]):
                sim[0] = t_course
            elif label == str(labels[1]):
                sim[1] = t_course
            elif label == str(labels[2]):
                sim[2] = t_course

        exp =  [[0] * 61] * 3
        for (columnName, columnData) in dfExp.items():
            label = list(columnData.iloc[0])
            
            time = np.linspace(0, 240, 61)
            if label == labels[0]:
                t_course = columnData.iloc[1:]
                exp[0] = [i/maxVal for i in t_course]
            elif label == labels[1]:
                t_course = columnData.iloc[1:]
                exp[1] = [i/maxVal for i in t_course]
            elif label == labels[2]:
                t_course = columnData.iloc[1:]
                exp[2] = [i/maxVal for i in t_course]
                
        err =  [[0] * 61] * 3
        for (columnName, columnData) in dfErr.items():
            label = list(columnData.iloc[0])
            
            time = np.linspace(0, 240, 61)
            if label == labels[0] and label != [5.0, 2.5, 0.001, 0, 90]:
                t_course = columnData.iloc[1:]
                err[0] = [i/maxVal for i in t_course]
            elif label == labels[1]:
                t_course = columnData.iloc[1:]
                err[1] = [i/maxVal for i in t_course]
            elif label == labels[2]:
                t_course = columnData.iloc[1:]
                err[2] = [i/maxVal for i in t_course]

                
        return sim, exp, err

    sim1, exp1, err1 = grabData(labels)

    if varyCondition == 'T7':
        varyIndex = 0
        vals = [1.0, 5.0, 20.0]
        colors = ['lightgreen', 'mediumseagreen', 'darkgreen']    
    elif varyCondition == 'RT':
        varyIndex = 1
        vals = [0.5, 2.5, 10.0]
        colors = ['lightsteelblue', 'royalblue', 'midnightblue']
    elif varyCondition == 'RNAse':
        varyIndex = 2
        vals = [0.001, 0.005, 0.02]
        colors = ['lightcoral', 'red', 'maroon']
        
    time = np.linspace(0, 240, 61)
 
    for i in range(0, len(exp1)):
        list_exp = exp1[i]
        list_err = err1[i]
        upper_y = []
        lower_y = []
        for j, val in enumerate(list_exp):
            upper_y.append(val + list_err[j])
            lower_y.append(val - list_err[j])

        ax1.fill_between(time, lower_y, upper_y, alpha = .2, color = colors[i])
        ax1.plot(time, exp1[i],  marker = None, linestyle = 'solid', color = colors[i])
        ax1.set_xscale('linear')
        
    for i in range(0, len(sim1)):
        ax3.plot(time, sim1[i],  marker = None, linestyle = 'dashed', color = colors[i])
    
    ax1.set_xlabel('Time (min)')
    ax1.set_ylabel('Normalized exp output')
    ax1.set_title('vRNA = 1fM', fontsize = 10, fontweight = 'bold')
    ax3.set_title('vRNA = 1fM', fontsize = 10, fontweight = 'bold')
    ax4.set_title('vRNA = 10fM', fontsize = 10, fontweight = 'bold')

    for i in range(0, len(labels)):
        labels[i][3] = 10
    
    sim10, exp10, err10 = grabData(labels)
    
    for i in range(0, len(sim1)):
        if "rep1" in dataID:
            if varyCondition == 'RT' and i == 0: #condition dropped due to high error
                continue
            else:
                ax4.plot(time, sim10[i],  marker = None, linestyle = 'dashed', color = colors[i])
        else:
            ax4.plot(time, sim10[i],  marker = None, linestyle = 'dashed', color = colors[i])
        
    for i in range(0, len(exp10)):
        list_exp = exp10[i]
        list_err = err10[i]
        upper_y = []
        lower_y = []
        for j, val in enumerate(list_exp):
            upper_y.append(val + list_err[j])
            lower_y.append(val - list_err[j])

        ax2.fill_between(time, lower_y, upper_y, alpha = .2, color = colors[i])
        ax2.plot(time, exp10[i],  marker = None, linestyle = 'solid', color = colors[i])
        ax2.set_xscale('linear')
    
    ax2.set_xlabel('Time (min)')
   
    ax2.set_title('vRNA = 10fM', fontsize = 10, fontweight = 'bold')
  
    if varyCondition == 'RT':
        objective = 4
        y_max_1 = y_max_RT[0]
        y_max_10 = y_max_RT[1]

    elif varyCondition == 'T7':
        objective = 5
        y_max_1 = y_max_T7[0]
        y_max_10 = y_max_T7[1]
     
    elif varyCondition == 'RNAse':
        objective = 6
        y_max_1 = y_max_RNase[0]
        y_max_10 = y_max_RNase[1]
        
    ax1.set_ylim(0, y_max_1)
    ax2.set_ylim(0, y_max_10)
    ax3.set_ylim(0, y_max_1)
    ax4.set_ylim(0, y_max_10)
    plt.savefig('./modeling_objective_' + str(objective) + '.svg', dpi = 600, bbox_inches="tight")

def plotModelingObjectives456(
        df_sim: pd.DataFrame, df_data: pd.DataFrame,
        df_error: pd.DataFrame, dataID: str, maxVal: float,
        y_max_RT: list, y_max_T7: list, y_max_RNase: list
) -> None:
    """Plots selected readout time courses for objectives 4, 5, 
        and 6 (enzyme sweeps). Calls resultsPanel for each 
        objective to generate plots
    
    Parameters
    ---------- 
    dfSim
        a dataframe containing the simulated data
   
    Returns
    -------
    None

    """
    cas = 90
    labels = [[1.0, 2.5, 0.005, 1, cas], [5.0, 2.5, 0.005, 1, cas], [20.0, 2.5, 0.005, 1, cas]]
    varyCondition = 'T7'
    resultsPanel(
        df_sim, df_data, df_error, dataID, labels, varyCondition,
        maxVal, y_max_RT, y_max_T7, y_max_RNase
    )

    labels = [[5.0, 0.5, 0.005, 1, cas], [5.0, 2.5, 0.005, 1, cas], [5.0, 10.0, 0.005, 1, cas]]
    varyCondition = 'RT'
    resultsPanel(
        df_sim, df_data, df_error, dataID, labels, varyCondition,
        maxVal, y_max_RT, y_max_T7, y_max_RNase
    ) 
    
    labels = [[5.0, 2.5, 0.001, 1, cas], [5.0, 2.5, 0.005,1, cas], [5.0, 2.5, 0.02, 1, cas]]
    varyCondition = 'RNAse'
    resultsPanel(
        df_sim, df_data, df_error, dataID, labels, varyCondition,
        maxVal, y_max_RT, y_max_T7, y_max_RNase
    )

def parityPlot(sim: list, exp: list) -> None:
    
    """
    Plots the experimental and simulated data in the form of a parity plot
    
    Parameters
    ---------- 
    sim
        a list of floats containing simulated values

    exp
        a list of floats containing experimental values
        
    Returns
    None
        
    Figures

    'fit_to_training_data_parity_plot.svg'
        a parity plot of the experimental vs simulated data
        
    """
    color_ = 'black'
        
    fig = plt.figure(figsize = (3.375,3))
    plt.plot(sim, exp, marker = 'o', markersize = 1, linestyle = 'none', color = color_)
    plt.ylabel('Experimental value')
    plt.xlabel('Simulated value')
    x_ = [0, 1]
    y_ = [0, 1]
    plt.plot(x_, y_, linestyle = '-', marker = None, color = 'dimgrey')
    plt.xticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    plt.yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    plt.savefig('./fit_to_training_data_parity_plot.svg', dpi = 600, bbox_inches="tight")