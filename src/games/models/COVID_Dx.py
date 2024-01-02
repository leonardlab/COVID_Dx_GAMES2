#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 09:11:29 2022

@author: kate
"""

import math
from typing import Tuple, List
import numpy as np
from scipy.integrate import odeint
from games.plots.plots_training_data import plot_training_data_2d
from games.models.Model_ODE_solver import ODE_solver, ODE_solver_D


class COVID_Dx:
    """
    Representation of COVID_Dx model

    """

    def __init__(
        self,
        parameters: List[float] = [1, 1, 1, 1, 1, 1],
        inputs: List[float] = None,
        mechanismID: str = "default",
        ode_solver_tolerance: str = "low"
    ) -> None:

        """Initializes COVID-Dx model.

        Parameters
        ----------
        parameters
            List of floats defining the parameters

        inputs
            List of floats defining the inputs

        mechanismID
            a string defining the mechanism identity

        Returns
        -------
        None

        """
        self.state_labels = [
            "vRNA (input)",
            "ssDNA p1",
            "ssDNA p2",
            "ssDNA p1:vRNA",
            "ssDNA p2:tRNA",
            "ssDNA p1:cvRNA",
            "ssDNA p2:ctRNA",
            "RT",
            "RNase H",
            "RT-ssDNA p1:vRNA",
            "RT-ssDNA p2:tRNA",
            "RT-ssDNA p1:cvRNA",
            "RT-ssDNA p2:ctRNA",
            "cDNA1:vRNA",
            "cDNA2:tRNA",
            "cDNA1:vRNA: RNase H",
            "cDNA2:tRNA: RNase H",
            "cDNA1:vRNA frag",
            "cDNA2:tRNA frag",
            "cDNA1:ssDNA p2",
            "cDNA2:ssDNA p1",
            "cDNA1:ssDNA p2:RT",
            "cDNA2:ssDNA p1:RT",
            "T7 RNAP",
            "dsDNA T7 target",
            "T7: dsDNA T7 target",
            "tRNA (target)",
            "target iCas13a-gRNA",
            "target aCas13a-gRNA",
            "dsRNA (input:target)",
            "quench-ssRNA-fluoro",
            "quencher",
            "fluorophore (output)"
        ]
        self.parameters = parameters
        self.inputs = inputs
        self.mechanismID = mechanismID
        self.ode_solver_tolerance = ode_solver_tolerance
        number_of_states = len(self.state_labels)
        x_init = np.zeros(number_of_states)
        self.initial_conditions = x_init

    def solve_single(self
    )-> Tuple[np.ndarray, list, np.ndarray]:
        """Solves COVID_Dx model for a single set of parameters and 
        inputs (doses)

        Parameters
        ----------
        none

        Returns
        -------
        solution
            An array of ODE solutions (rows are timepoints and columns are model states)

        t
            A 1D array of time values corresponding to the rows in solution

        """

        if self.mechanismID == "D":
            solver = ODE_solver_D()

            C_scale = 10 ** 6

            self.initial_conditions[0] = self.inputs[3] * .000001  # x_v
            self.initial_conditions[0] = self.initial_conditions[0] * C_scale #x_v' sent into ODEs (x_v' = x_v * 10^6)   
            self.initial_conditions[1] = 250 # x_p1
            self.initial_conditions[2] = 250 # x_p2
            self.initial_conditions[7] = self.inputs[1] * 139.1 # x_RT
            self.initial_conditions[8] = self.inputs[2] * 6060 # x_RNase
            self.initial_conditions[23] = self.inputs[0] * 16.16 # x_T7
            self.initial_conditions[27] = self.inputs[4]/2 # x_iCas13
            self.initial_conditions[30] = 2500 # x_qRf
            solver.set_initial_condition(np.array(self.initial_conditions))

            #Parameters
            k_cas13  = self.parameters[0] #nM-1 min-1
            k_degv = self.parameters[1] #nM-1 min-1
            k_txn = self.parameters[2] #min-1
            k_FSS = self.parameters[3] #min-1
            a_RHA = self.parameters[4]
            b_RHA = self.parameters[5]
            c_RHA = self.parameters[6]

            k_bds = k_cas13 #nM-1 min-1
            k_RTon = 0.024 #nM-1 min-1
            k_RToff = 2.4 #min-1
            k_T7on = 3.36 #nM-1 min-1
            k_T7off = 12 #min-1
            k_SSS = k_FSS #min-1
            k_degRrep = k_degv  #nM-1 min-1
            k_RNaseon = 0.024 #nM-1 min-1
            k_RNaseoff = 2.4 #min-1
            the_rates = np.array([k_degv, k_bds, k_RTon, k_RToff, k_RNaseon, k_RNaseoff, 
                                k_T7on, k_T7off, k_FSS, a_RHA, b_RHA, c_RHA, k_SSS, 
                                k_txn, k_cas13, k_degRrep]).astype(float)
            solver.set_rates(the_rates)
            solver.abs_tol = 1e-13
            solver.rel_tol = 1e-10
            solver.complete_output = 0
            solver.conservation_form = True
            solver.dist_type = "expon"
        
            #Time-stepping
            timesteps = (240 * 100) + 1
            final_time = 240
            tspace = np.linspace(0, final_time, timesteps)
            
            #Set solver algorithm
            solver.solver_alg = "LSODA"
            solver.ode_solver_tolerance = self.ode_solver_tolerance
            solver.k_loc_deactivation = self.parameters[7]
            solver.k_scale_deactivation = self.parameters[8]

            solver.mechanism_B = "yes"
            solver.mechanism_C = "yes"
            solver.txn_poisoning = "no"
            
        else:
            solver = ODE_solver()
            
            C_scale = 10 ** 6
            self.initial_conditions[0] = self.inputs[3] * .000001  # x_v
            self.initial_conditions[0] = self.initial_conditions[0] * C_scale #x_v' sent into ODEs (x_v' = x_v * 10^6)   
            self.initial_conditions[1] = 250 # x_p1
            self.initial_conditions[2] = 250 # x_p2
            self.initial_conditions[7] = self.inputs[1] * 139.1 # x_RT
            self.initial_conditions[8] = self.inputs[2] * 6060 # x_RNase
            self.initial_conditions[23] = self.inputs[0] * 16.16 # x_T7
            self.initial_conditions[27] = self.inputs[4]/2 # x_iCas13
            self.initial_conditions[30] = 2500 # x_qRf
            solver.set_initial_condition(np.array(self.initial_conditions))

            #Parameters
            k_cas13  = self.parameters[0] #nM-1 min-1
            k_degv = self.parameters[1] #nM-1 min-1
            k_txn = self.parameters[2] #min-1
            k_FSS = self.parameters[3] #min-1
            k_RHA = self.parameters[4] #min-1
            
            k_bds = k_cas13 #nM-1 min-1
            k_RTon = .024 #nM-1 min-1
            k_RToff = 2.4 #min-1
            k_T7on = 3.36 #nM-1 min-1
            k_T7off = 12 #min-1
            k_SSS = k_FSS #min-1
            k_degRrep = k_degv  #nM-1 min-1
            k_RNaseon = .024 #nM-1 min-1
            k_RNaseoff = 2.4 #min-1
            the_rates = np.array([k_degv, k_bds, k_RTon, k_RToff, k_RNaseon, k_RNaseoff,
                                k_T7on, k_T7off, k_FSS, k_RHA, k_SSS, k_txn, k_cas13,
                                k_degRrep]).astype(float)
            solver.set_rates(the_rates)
            solver.abs_tol = 1e-13
            solver.rel_tol = 1e-10
            solver.complete_output = 0
            solver.conservation_form = True
            solver.dist_type = "expon"

            #Time-stepping
            timesteps = (240 * 100) + 1
            final_time = 240
            tspace = np.linspace(0, final_time, timesteps)
            
            #Set solver and algorithm
            solver.solver_alg = "LSODA"
            solver.ode_solver_tolerance = self.ode_solver_tolerance
            solver.k_loc_deactivation = self.parameters[5]
            solver.k_scale_deactivation = self.parameters[6]
            
            
            if self.mechanismID == "A":
                solver.mechanism_B = "no"
                solver.mechanism_C = "no"
                solver.txn_poisoning = "no"
            
            elif self.mechanismID == "B":
                solver.mechanism_B = "yes"
                solver.mechanism_C = "no"
                solver.txn_poisoning = "no"
            
            elif self.mechanismID == "C":
                solver.mechanism_B = "yes"
                solver.mechanism_C = "yes"
                solver.txn_poisoning = "no"

        #Solve equations
        solution, t = solver.solve(tspace)
        
        #Round results
        solution = np.around(solution, decimals = 10)
        
        #Define the time course of the readout
        timecourse_readout = solution[:, -1]   
    
        #Unscale vRNA to original units
        vRNA_unscaled = [i/C_scale for i in solution[:,0]] #vRNA = vRNA' / 1000000
        solution[:,0] = vRNA_unscaled
        
        #Restructure t and readout time course to match exp data sampling
        t = t[::400]
        timecourse_readout = timecourse_readout[::400]

        #Restructure other model state time courses to match exp data and t
        solutions = [solution[:, i][::400] for i in range(0, np.shape(solution)[-1])]

        return t, solutions, timecourse_readout

    def solve_experiment(self, x: list) -> list[float]:
        """Solve synTF_Chem model for a list of ligand values.

        Parameters
        ----------
        x
            a list of floats containing the independent variable

        dataID
            a string defining the dataID

        parameter_labels
            a list of strings defining the parameter labels

        Returns
        -------
        solutions
            A list of floats containing the value of the reporter protein
            at the final timepoint for each ligand amount

        """

        solutions = []

        for doses in x:
            self.inputs = doses
            t, _, reporter_timecourse = self.solve_single()

            if len(reporter_timecourse) == len(t):
                reporter_timecourse = reporter_timecourse
            
            else:
                reporter_timecourse = [0] * len(t)
            
            for i in reporter_timecourse:
                solutions.append(float(i))

        return solutions

    @staticmethod
    def normalize_data(solutions_raw: List[float]) -> List[float]:
        """Normalizes data by maximum value

        Parameters
        ----------
        solutions_raw
            a list of floats defining the solutions before normalization

        dataID
            a string defining the dataID

        Returns
        -------
        solutions_norm
            a list of floats defining the dependent variable for the given
            dataset (after normalization)

        """

        if max(solutions_raw) == 0:
            solutions_norm = [0]*len(solutions_raw)

        else:
            solutions_norm = [i/max(solutions_raw) for i in solutions_raw]

        return solutions_norm

    # @staticmethod
    # def plot_training_data(
    #     x: list,
    #     solutions_norm: List[float],
    #     exp_data: List[float],
    #     exp_error: List[float],
    #     filename: str,
    #     run_type: str,
    #     context: str,
    #     dataID: str,
    # ) -> None:
    #     """
    #     Plots training data and simulated training data for a single parameter set

    #     Parameters
    #     ----------
    #     x
    #         list of floats defining the independent variable

    #     solutions_norm
    #         list of floats defining the simulated dependent variable

    #     exp_data
    #         list of floats defining the experimental dependent variable

    #     exp_error
    #         list of floats defining the experimental error for the dependent variable

    #     filename
    #        a string defining the filename used to save the plot

    #     run_type
    #         a string containing the data type ('PEM evaluation' or else)

    #     context
    #         a string defining the file structure context

    #     dataID
    #         a string defining the dataID

    #     Returns
    #     -------
    #     None"""

    #     # define plot settings
    #     if run_type == "default":
    #         plot_color = "black"
    #         marker_type = "o"

    #     elif run_type == "PEM evaluation":
    #         plot_color = "dimgrey"
    #         marker_type = "^"

    #     if dataID == "ligand dose response":
    #         y_label = "Rep. protein (au)"
    #         x_label = "Ligand (nM)"
    #         x_scale = "symlog"
    #         plot_settings = x_label, y_label, x_scale, plot_color, marker_type
    #         plot_training_data_2d(
    #             x, solutions_norm, exp_data, exp_error, filename, plot_settings, context
    #         )

    #     elif dataID == "ligand dose response and DBD dose response":
    #         # Define plot settings for ligand dose response
    #         y_label = "Rep. protein (au)"
    #         x_label = "Ligand (nM)"
    #         x_scale = "symlog"
    #         plot_settings = x_label, y_label, x_scale, plot_color, marker_type

    #         # plot ligand dose response
    #         filename_1 = filename + "ligand dose response"
    #         plot_training_data_2d(
    #             x[:11],
    #             solutions_norm[:11],
    #             exp_data[:11],
    #             exp_error[:11],
    #             filename_1,
    #             plot_settings,
    #             context,
    #         )

    #         # Define plot settings for DBD dose response
    #         y_label = "Rep. protein (au)"
    #         x_label = "DBD plasmid dose (ng)"
    #         x_scale = "linear"
    #         plot_settings = x_label, y_label, x_scale, plot_color, marker_type

    #         # Plot DBD dose response @ 20ng AD
    #         filename_2 = filename + "DBD dose response 20ng AD"
    #         plot_training_data_2d(
    #             x[11:19],
    #             solutions_norm[11:19],
    #             exp_data[11:19],
    #             exp_error[11:19],
    #             filename_2,
    #             plot_settings,
    #             context,
    #         )

    #         # Plot DBD dose response @ 10ng AD
    #         filename_3 = filename + "DBD dose response 10ng AD"
    #         plot_training_data_2d(
    #             x[19:],
    #             solutions_norm[19:],
    #             exp_data[19:],
    #             exp_error[19:],
    #             filename_3,
    #             plot_settings,
    #             context,
    #         )
