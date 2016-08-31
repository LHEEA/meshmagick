#!/usr/bin/env python
#  -*- coding: utf-8 -*-
"""This module contains the base classes of the Framework

Base classes are grouped in a central module so as to avoid some import problems
and circular dependencies

"""
from itertools import count

# TODO: Ajouter une transformation pour toutes les frames ! Pour Ground et Null, pas de transformation,
# ou alors l'identite...

# -------------------------------------------------------
# Design Patterns
# -------------------------------------------------------
class _Singleton(type):
    _instances = {}

    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            cls._instances[cls] = super(_Singleton, cls).__call__(*args, **kwargs)
        return cls._instances[cls]

# -------------------------------------------------------
# Frames
# -------------------------------------------------------
_frames_ids = count(0)

class _BaseFrame(object):
    """The base class to represent a Frame in 3D space

    """

    def __init__(self, name=None, parent_frame='ground'):
        """

        Parameters
        ----------
        name : str
            The name of the frame object
        parent_frame : _BaseFrame
            The parent Frame of the current frame

        """
        self._id = _frames_ids.next()
        # TODO: travailler les conventions de nommage pour un nom plus explicite
        if not name:
            name = '_'.join(('Frame', str(self._id)))

        self._name = str(name)
        self.parent_frame = parent_frame

    def __str__(self):
        str_repr = "Frame '%s'" % self._name
        return str_repr

    __repr__ = __str__

    @property
    def name(self):
        """
        Return the name of the Frame

        Returns
        -------
        name : str
            The name of the Frame
        """
        return self._name

    @name.setter
    def name(self, _name):
        self._name = str(_name)

    @property
    def parent_frame(self):
        """
        Return the parent Frame of the current frame

        Returns
        -------
        parent_frame : _BaseFrame
            The Frame's parent frame
        """
        return self._parent_frame

    @parent_frame.setter
    def parent_frame(self, _parent_frame):
        if not isinstance(_parent_frame, (_BaseFrame, _BodyReferenceFrame, _GroundFrame)):
            raise TypeError("parent_frame of BaseFrame must be a BaseFrame object")
        self._parent_frame = _parent_frame
        return


# TODO: implementer une methode __null__ pour permettre de faire 'while Frame:'
class _NullFrame(object):
    __metaclass__ = _Singleton

    # def __init__(self):
    #     pass

    def __str__(self):
        return "Frame 'Null_Frame'"

    __repr__ = __str__

    @property
    def name(self):
        return 'Null_Frame'

    @property
    def parent_frame(self):
        return None


null_frame = _NullFrame()


# FIXME: Singleton necessaire ???
class _GroundFrame(_BaseFrame):
    __metaclass__ = _Singleton

    def __init__(self):
        self._id = _frames_ids.next()

    def __str__(self):
        return "Frame 'Ground_Frame'"

    __repr__ = __str__

    @property
    def name(self):
        return 'Ground_Frame'

    @property
    def parent_frame(self):
        return null_frame


ground_frame = _GroundFrame()


# TODO: Est-il necessaire de definir une classe BodyFixedFrame ? Qu'est-ce que cela apporte par rapport a BaseFrame ?
#  --> ca ajoute body_owner !!!!

class BodyFixedFrame(_BaseFrame):
    """Defines a frame that is attached (belongs) to a body

    """

    def __init__(self, body_owner, name=None):
        if not isinstance(body_owner, _BaseBody):
            raise TypeError('body owner of a BodyFixedFrame must be a _BaseBody object')
        self._body_owner = body_owner

        super(BodyFixedFrame, self).__init__(name=name, parent_frame=body_owner.reference_frame)


    def __str__(self):
        str_repr = "Frame '%s' attached to body %s" % (self.name, self._body_owner.name)
        return str_repr

    __repr__ = __str__

    @property
    def body_owner(self):
        return self._body_owner


# TODO: BodyReferenceFrame devrait deriver de BodyFixedFrame, ca ajoute le champ children ...
class _BodyReferenceFrame(BodyFixedFrame):
    def __init__(self, body_owner, parent_frame=ground_frame):
        super(BodyFixedFrame, self).__init__(name=None, parent_frame=parent_frame)

        # self._id = _frames_ids.next()
        # self._body_owner = body_owner
        # self._parent_frame = ground_frame
        self._children = dict()

    def __str__(self):
        return "Reference frame of body '%s'" % self._body_owner.name

    __repr__ = __str__

    # We redefine props name and body_owner so as to prevent name modification
    # from the outside
    @property
    def name(self):
        return '_'.join((self._body_owner.name, 'Reference_Frame'))

    @property
    def body_owner(self):
        return self._body_owner

    @property
    def parent_frame(self):
        return self._parent_frame

    @parent_frame.setter
    def parent_frame(self, _parent_frame):
        if not isinstance(_parent_frame, (_BaseFrame, _GroundFrame)):
            raise TypeError, "Parent Frame of a BodyReferenceFrame must be a _BaseFrame or a _GroundFrame object"
        self._parent_frame = _parent_frame

    @property
    def children(self):
        return self._children

    def add_child(self, frame):
        frame_name = frame.name
        if self._children.has_key(frame_name):
            err_msg = "Naming clash : A frame with name %s has already " \
                      "been added referenced with respect to %s" % (frame_name, self.name)
            raise ValueError(err_msg)
        self._children[frame_name] = frame
        return


# -------------------------------------------------------
# Bodies
# -------------------------------------------------------
_bodies_ids = count(0)

# TODO: il faut un comptage local des frames !!!
class _BaseBody(object):

    def __init__(self, name=None):
        """

        Parameters
        ----------
        name : str, optional
            The body's name

        """
        self._id = _bodies_ids.next()
        if not name:
            name = '_'.join(('Body', str(self._id)))
        self._name = str(name)

        self._reference_frame = _BodyReferenceFrame(self)


    def __str__(self):
        str_repr = "Body '%s'" % self._name
        return str_repr

    __repr__ = __str__

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, _name):
        self._name = str(_name)
        return

    @property # No setter as we want to prevent reference frame change from the outside
    def reference_frame(self):
        return self._reference_frame

    def add_fixed_frame(self, name=None):

        new_frame = BodyFixedFrame(self, name=name)
        self._reference_frame.add_child(new_frame)

        return new_frame

    @property
    def frames(self):
        frames = self.reference_frame.children.copy()
        frames['reference_frame'] = self.reference_frame
        return frames

# FIXME: si tous les arguments sont statiques, pourquoi utiliser le pattern singleton ? --> inutile
class _NullBody(object):
    """
    A body class representing an undefined body
    """
    __metaclass__ = _Singleton

    # def __init__(self):
    #     pass

    def __str__(self):
        return "Body 'Null_Body'"

    __repr__ = __str__

    @property
    def name(self):
        return 'Null_Body'

    @property
    def reference_frame(self):
        return None


_null_body = _NullBody()


class _GroundBody(object):
    __metaclass__ = _Singleton

    def __init__(self):
        self._id = _bodies_ids.next()

    def __str__(self):
        return "Body 'Ground_Body'"

    @property
    def name(self):
        return 'Ground_Body'

    @property
    def reference_frame(self):
        return ground_frame


ground_body = _GroundBody()
