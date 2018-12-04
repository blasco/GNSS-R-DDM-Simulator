#!/usr/bin/env python3

from gnssr.tds.tds_data import *
from gnssr.targets import *

# Di Simone Oil Platform Data
file_root_name = 'raw/L1B/2015-04-01-H00'
target = targets['hibernia']
group = '000095'
index_start = 525
tds = tds_data(file_root_name)

import socket, pickle
import time

HOST = '127.0.0.1'  # The server's hostname or IP address
PORT = 65432        # The port used by the server

for i in range(index_start-200, index_start+100):
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.connect((HOST, PORT))
        ddm = tds.rootgrp.groups[group].variables['DDM'][i].data
        s.sendall(pickle.dumps(ddm))
        s.close()
    if i < 525:
        time.sleep(0.1)
    else:
        time.sleep(1)
