#!/usr/bin/env python3

import struct
import socket
import pickle
import sys
import traceback
from threading import Thread

from gnssr.tds.tds_data import *
from gnssr.targets import *

import tm

tm_processor = tm.processor()

def send_msg(sock, msg):
    # Prefix each message with a 4-byte length (network byte order)
    msg = struct.pack('>I', len(msg)) + msg
    sock.sendall(msg)

def main():
    start_server()

def start_server():
    host = "127.0.0.1"
    port = 8888         # arbitrary non-privileged port

    soc = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    soc.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)   # SO_REUSEADDR flag tells the kernel to reuse a local socket in TIME_WAIT state, without waiting for its natural timeout to expire
    print("Socket created")

    try:
        soc.bind((host, port))
    except:
        print("Bind failed. Error : " + str(sys.exc_info()))
        sys.exit()

    soc.listen(5)      # queue up to 5 requests
    print("Socket now listening")

    # infinite loop- do not reset for every requests
    while True:
        connection, address = soc.accept()
        ip, port = str(address[0]), str(address[1])
        print("Connected with " + ip + ":" + port)

        try:
            Thread(target=client_thread, args=(connection, ip, port)).start()
        except:
            print("Thread did not start.")
            traceback.print_exc()

    soc.close()

def client_thread(connection, ip, port, max_buffer_size = 5120):
    is_active = True

    while is_active:
        client_input = receive_input(connection, max_buffer_size)

        if client_input == "--QUIT--":
            print("Client is requesting to quit")
            connection.close()
            print("Connection " + ip + ":" + port + " closed")
            is_active = False

        elif client_input == 'GET_TM':
            index = tm_processor.find_new_data()
            if index > 0:
                connection.sendall("NEW_DATA".encode("utf8"))
                filename = os.path.join(os.environ['PROJECT_SRC_ROOT'],'test_setup/obc/fs/ddm___{0}'.format(index))
                with open(filename, 'rt') as f:
                    data = []
                    for line in f:
                        data.append(line)
                    send_msg(connection, pickle.dumps(data))
            else:
                connection.sendall("NO_NEW_DATA".encode("utf8"))

        else:
            print("Unknown command: {0}".format(client_input))
            connection.sendall("-".encode("utf8"))

def receive_input(connection, max_buffer_size):
    client_input = connection.recv(max_buffer_size)
    client_input_size = sys.getsizeof(client_input)

    if client_input_size > max_buffer_size:
        print("The input size is greater than expected {}".format(client_input_size))

    decoded_input = client_input.decode("utf8").rstrip()  # decode and strip end of line
    return decoded_input

if __name__ == "__main__":
    main()
