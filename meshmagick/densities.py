#!/usr/bin/env python
#  -*- coding: utf-8 -*-
"""Module to retrieve some densities of different mediums"""


_DENSITIES = {'CONCRETE': 2300.,
              'REINFORCED_CONCRETE': 2400.,
              'FRESH_WATER': 1000.,
              'SALT_WATER': 1025.,
              'SAND': 1600.,
              'STEEL': 7850.,
              'ALUMINUM': 2700.,
              'LEAD': 11350.,
              'TITANIUM': 4500.,
              'POLYSTYRENE': 1050.,
              'GASOLINE': 750.,
              'DIESEL_FUEL': 850.,
              'ETHANOL': 789.,
              'AIR_20DEGC': 1.204,
              'BUTANE': 2.7,
              'PROPANE': 2.01,
              'HYDROGEN_-252DEGC': 70.,
              'NITROGEN_-195DEGC': 810.
              }


def get_density(medium):
    try:
        density = _DENSITIES[medium]
    except KeyError:
        raise KeyError('Medium %s not known...' % medium)
    return density

def list_mediums():
    return _DENSITIES.keys()
