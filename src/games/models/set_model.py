#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  3 14:57:05 2022

@author: kate
"""
import json
from games.models.COVID_Dx import COVID_Dx



def set_model():
    """
    Defines the model class to use depending on the modelID defined in Settings

    Parameters
    ----------
    None

    Returns
    -------
    model
        object defining the model

    """

    given_model = COVID_Dx(
        parameters=settings["parameters"],
        mechanismID=settings["mechanismID"]
    )


    return given_model


file = open("/Users/kdreyer/Documents/Github/COVID_Dx_GAMES2/src/games/config/config_COVID_Dx_D.json", encoding="utf-8")
settings = json.load(file)
model = set_model()
