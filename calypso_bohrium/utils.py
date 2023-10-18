#!/usr/bin/env python

import os


def get_value(key):
    # temp = os.popen(f"grep {key} input.dat -i ").read().split("=")[-1].strip()
    temp = os.popen(f"grep {key} input.dat -i ").read().split("#")[0].split("=")[-1].strip()
    return temp

