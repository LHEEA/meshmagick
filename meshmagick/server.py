#!/usr/bin/env python
#  -*- coding: utf-8 -*-

from multiprocessing import Process
from time import sleep
import zmq
import numpy as np
import sys
import Queue


def send_array(socket, A, flags=0, copy=True, track=False):
    """send a numpy array with metadata"""
    md = dict(
        dtype = str(A.dtype),
        shape = A.shape,
    )
    socket.send_json(md, flags|zmq.SNDMORE)
    return socket.send(A, flags, copy=copy, track=track)

def recv_array(socket, flags=0, copy=True, track=False):
    """recv a numpy array"""
    md = socket.recv_json(flags=flags)
    msg = socket.recv(flags=flags, copy=copy, track=track)
    buf = buffer(msg)
    A = np.frombuffer(buf, dtype=md['dtype'])
    return A.reshape(md['shape'])





def client(port):
    context = zmq.Context()
    print "Connecting to server with port %s" % port
    socket = context.socket(zmq.PULL)
    socket.connect("tcp://localhost:%u" % port)
    
    msgs = Queue.Queue()
    
    while True:
        # Cette boucle sera celle du iren.Start() dans un callback...
        # arr = recv_array(socket)
        msg = socket.recv_string()
        msgs.put(msg)
        print 'client ', msg
        
        if msg == 'fin':
            break
    
    while not msgs.empty():
        msg = msgs.get()
        print msg
        


if __name__ == "__main__":
    
    # Running server
    port = 5561
    
    context = zmq.Context()
    socket = context.socket(zmq.PUSH)
    socket.bind("tcp://*:%u" % port)
    print "Running server on port: ", port
    
    
    print "Starting client"
    Process(target=client, args=(port,)).start()
    
    txt = 'Ma pauvre choupi est toute malade en cette fin de journee'.split()
    
    for msg in txt:
        try:
            # send_array(socket, arr, copy=True)
            socket.send_string(msg)
        except:
            print 'fail to send'
        sleep(0.1)
    # Now we can connect a client to all these servers
    
