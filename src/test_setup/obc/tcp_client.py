#!/usr/bin/env python3

import struct
import socket
import pickle
import sys

import target

def recvall(sock, n):
    # Helper function to recv n bytes or return None if EOF is hit
    data = b''
    while len(data) < n:
        packet = sock.recv(n - len(data))
        if not packet:
            return None
        data += packet
    return data

def recv_msg(sock):
    # Read message length and unpack it into an integer
    raw_msglen = recvall(sock, 4)
    if not raw_msglen:
        return None
    msglen = struct.unpack('>I', raw_msglen)[0]
    # Read the message data
    return recvall(sock, msglen)

def main():
    target_processor = target.processor()

    soc = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    host = "127.0.0.1"
    port = 8889

    try:
        soc.connect((host, port))
    except:
        print("Connection error")
        sys.exit()

    print("Enter 'quit' to exit")
    command = input(" -> ")

    while command != 'quit':

        if command == "START":
            for index in range(0, 97):
                print(index)
                soc.sendall("GET_DDM".encode("utf8"))
                data = recv_msg(soc)
                ddm = pickle.loads(data)

                soc.sendall("GET_METADATA".encode("utf8"))
                data = recv_msg(soc)
                metadata = pickle.loads(data)

                soc.sendall("NEXT_DDM".encode("utf8"))

                target_processor.load_ddm(ddm, metadata)
                target_processor.process_ddm(index)

                #target_processor.plot_targets()
        else:
            soc.sendall(command.encode("utf8"))
            res = soc.recv(5120).decode("utf8")
            if res == "-":
                pass        # null operation

        command = input(" -> ")

    soc.send(b'--QUIT--')

if __name__ == "__main__":
    main()
