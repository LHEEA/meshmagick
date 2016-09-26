#!/usr/bin/env python
#  -*- coding: utf-8 -*-
import zmq
from time import sleep



class Client():
    def __init__(self):
        print 'client instantiation'
        
    def __call__(self, *args, **kwargs):
        print 'client called'

