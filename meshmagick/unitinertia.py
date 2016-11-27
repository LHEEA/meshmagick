#!/usr/bin/env python
#  -*- coding: utf-8 -*-

import numpy as np

# FIXME: attention, changer les signes pour les produits d'inertie !

# TODO: ajouter la production d'inerties de solides connus --> utile pour comparaison !!


class UnitInertia(object):
    def __init__(self, volume, cog, Ixx, Iyy, Izz, Iyz, Ixz, Ixy, point=None):
        
        self._volume = float(volume)
        
        assert len(cog) == 3
        self._cog = np.asarray(cog, dtype=np.float)
        
        if point is None:
            self._point = self._cog
        else:
            assert len(point) == 3
            self._point = np.asarray(point, dtype=np.float)
        
        # We only store coefficients
        self._inertia_coeffs = [Ixx, Iyy, Izz, Iyz, Ixz, Ixy]
        
        self._corr_G = self._huygens_correction(self._point)
        

    @property
    def is_at_cog(self):
        return np.all(self._point == self._cog)
    
    def get_inertia_mat(self, density=1., at_cog=False):
        
        i = self._inertia_coeffs
        inertia_mat = np.array([[i[0], -i[5], -i[4]],
                                [-i[5], i[1], -i[3]],
                                [-i[4], -i[3], i[2]]], dtype=np.float)
        
        if at_cog and not self.is_at_cog:
            inertia_mat -= self._corr_G
        
        return inertia_mat * density
    
    @property
    def reduction_point(self):
        return self._point
    
    @property
    def cog(self):
        return self._cog
    
    @reduction_point.setter
    def reduction_point(self, point):
        """Change the reduction point"""
        assert len(point) == 3
        point = np.asarray(point, dtype=np.float)
        inertia_mat = self.get_inertia_mat(at_cog=True)
        
        self._point = point
        self._corr_G = self._huygens_correction(point)
        
        inertia_mat += self._corr_G
        self._inertia_coeffs = inertia_mat.reshape((9,))[[0, 4, 8, 5, 2, 1]] * [1, 1, 1, -1, -1, -1]
        return
    
    def set_point_at_cog(self):
        self.reduction_point = self._cog
        return
    
    def _huygens_correction(self, point):
        p_g = self._cog - point
        return self._volume * (np.dot(p_g, p_g)*np.eye(3) - np.outer(p_g, p_g))
    
    def get_mass(self, density=1.):
        return self._volume * density
    
    
    def get_coeffs(self, density=1.):
        
        cl = self._inertia_coeffs * density
        coeffs = {'Ixx': cl[0],
                  'Iyy': cl[1],
                  'Izz': cl[2],
                  'Iyz': cl[3],
                  'Ixz': cl[4],
                  'Ixy': cl[5],
                  }
        return coeffs


if __name__ == '__main__':
    inertia = UnitInertia(1, [0, 0, 0], 1, 2, 3, 4, 5, 6, point=[1, 2, 3])
    print inertia.get_inertia_mat()
    print inertia.reduction_point
    print inertia.get_inertia_mat(at_cog=True)
    
    inertia.reduction_point = [2, 2, 2]
    
    print inertia.get_inertia_mat()
    
    print inertia.get_inertia_mat(at_cog=True)
    
    inertia.set_point_at_cog()
    
    print inertia.get_inertia_mat()
    
    print inertia.get_inertia_mat(at_cog=True)
    
