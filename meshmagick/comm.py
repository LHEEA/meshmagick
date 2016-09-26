#!/usr/bin/env python
#  -*- coding: utf-8 -*-

from multiprocessing import Process
from time import sleep
import zmq
# from zmq.devices import ProcessDevice

def data_server(port):
    context = zmq.Context()
    socket = context.socket(zmq.SUB)
    socket.connect("tcp://127.0.0.1:%u" % port)
    
    print 'Data server started'
    topic_filter = 'topic'
    socket.setsockopt(zmq.SUBSCRIBE, topic_filter)
    while True:
        msg = socket.recv()
        print 'received', msg
        sleep(0.1)

    
    

if __name__ == '__main__':
    
    port = 6000
    
    context = zmq.Context()
    socket = context.socket(zmq.PUB)
    socket.setsockopt(zmq.SNDHWM, 10000)
    socket.bind("tcp://127.0.0.1:%u" % port)
    
    
    Process(target=data_server, args=(port,)).start()
    # Process(target=data_server, args=(port,)).start()
    # Process(target=data_server, args=(port,)).start()
    # Process(target=data_server, args=(port,)).start()
    #
    
    iter = 0
    while True:
        topic = 'topic'
        msg = '%s %u' % (topic, iter)
        print 'Sending ', msg
        socket.send(msg)
        sleep(1)
        
        iter += 1
