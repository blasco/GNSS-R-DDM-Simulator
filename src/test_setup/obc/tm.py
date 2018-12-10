#!/usr/bin/env python

import os

class processor:

    def __init__(self):
        self.last_sent_index = 0
        self.fs_path = os.path.join(os.environ['PROJECT_SRC_ROOT'],'test_setup/obc/fs')
        for filename in os.listdir(self.fs_path):
            os.remove(self.fs_path + '/' + filename)

    def find_new_data(self):
        for filename in os.listdir(self.fs_path):
            index = int(filename.split('___')[1])
            if index == self.last_sent_index + 1:
                self.last_sent_index = index
                return index 
            else:
                pass
        return 0

    def save_sent_index(self, index):
        self.last_sent_index = index
