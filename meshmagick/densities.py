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
    """Get the density of medium.
    
    Parameters
    ----------
    medium : str
        The medium name

    Returns
    -------
    float
        Density of medium (kg/m**3)
    """
    try:
        density = _DENSITIES[str(medium).upper()]
    except KeyError:
        raise KeyError('Medium %s not known...' % medium)
    return density


def list_medium():
    """Get the list of available medium.
    
    Returns
    -------
    list
        List of available medium
    """
    return _DENSITIES.keys()


# def get_table():
#
#     col_width = 22
#     hline = '+{0:s}+{0:s}+\n'.format('-' * col_width)
#     table = hline
#     table += '|{:<{n}s}|{:>{n}s}|\n'.format('NAME', 'DENSITY (KG/M**3)', n=col_width)
#     table += hline
#
#     for key in _DENSITIES:
#         table +=
#         table += hline
#
#     return table
#
# if __name__ == '__main__':
#
#     print get_table()
