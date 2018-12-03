#!/usr/bin/env python3

from gnssr.tds.tds_data import *
from gnssr.targets import *

file_root_name = 'raw/L1B/2017-11-13-H18'
target = targets['devils_tower']
group = '000057'
index = 426

tds = tds_data(file_root_name)
tds.set_group_index(group, index)

import socket, pickle

HOST = '127.0.0.1'  # The server's hostname or IP address
PORT = 65432        # The port used by the server

with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
    s.connect((HOST, PORT))
    ddm = tds.rootgrp.groups[tds.group].variables['DDM'][tds.index].data
    s.sendall(pickle.dumps(ddm))
    s.close()
