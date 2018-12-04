#!/usr/bin/env python

import socket, pickle

HOST = '127.0.0.1'  # Standard loopback interface address (localhost)
PORT = 65432        # Port to listen on (non-privileged ports are > 1023)

import matplotlib.pyplot as plt

from gnssr.tds.detection.find_targets import *

p = target_processor()
while True:
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.bind((HOST, PORT))
        s.listen()
        conn, addr = s.accept()
        data = []
        with conn:
            print('Connected by', addr)
            while True:
                packet = conn.recv(2000)
                if not packet: break
                data.append(packet)
            serialized_data = b''.join(data)
            ddm = pickle.loads(serialized_data)
            p.process_ddm(ddm)
