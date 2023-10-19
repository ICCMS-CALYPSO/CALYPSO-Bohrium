#!/usr/bin/env python

import os

def get_value(key):
    # temp = os.popen(f"grep {key} input.dat -i ").read().split("=")[-1].strip()
    temp = (
        os.popen(f"grep {key} input.dat -i ")
        .read()
        .split("\n")
    )
    for t in temp:
        if t.strip().startswith("#"):
            continue
        elif "=" in t:
            val = t.split("#")[0].split("=")[-1]
            break
    else:
        val = ""
    return val
