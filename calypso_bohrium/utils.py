#!/usr/bin/env python

import os


def get_value(key):
    temp = os.popen(f"grep {key} input.dat").read().split("=")[-1].strip()
    return temp

