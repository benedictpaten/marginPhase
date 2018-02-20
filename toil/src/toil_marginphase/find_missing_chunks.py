#!/usr/bin/env python
from __future__ import print_function
import argparse
import glob
import json
import sys
import numpy as np
import os

EXCLUDED_PARTS = ['merged']

def parse_args():
    parser = argparse.ArgumentParser("Consolidates output learned parameters for chunked marginPhase run")
    parser.add_argument('--input_directory', '-i', dest='input_directory', default=".", type=str,
                       help='What directory we\'re searching in')
    parser.add_argument('--depth', '-d', dest='depth', default=5, type=int,
                       help='Output file. Default:stdout')
    parser.add_argument('--find_chunks', '-f', dest='find_chunks', action='store_true', default=False,
                       help='Actually find missing chunks')
    parser.add_argument('--threshold', '-t', dest='threshold', default=4, type=str,
                       help='Min number of files to not count as failure')

    return parser.parse_args()


def add_file_to_name_tree(file_name_tree, file_parts):
    if file_parts is None or len(file_parts) == 0: return
    car = file_parts[0]
    cdr = file_parts[1:]
    if car not in file_name_tree:
        file_name_tree[car] = dict()
    add_file_to_name_tree(file_name_tree[car], cdr)


def get_leaf_count(file_name_tree):
    if file_name_tree is None or len(file_name_tree) == 0:
        return 1
    return sum([get_leaf_count(file_name_tree[name]) for name in file_name_tree.keys()])


def get_file_names_for_depth(file_name_tree, depth):
    if file_name_tree is None or len(file_name_tree) == 0:
        return [""]
    if depth == 0: return [".."]
    names = list(file_name_tree.keys())
    names.sort()

    output = list()
    for name in names:
        children = get_file_names_for_depth(file_name_tree[name], depth - 1)
        for child in children:
            output.append(name + "." + child if len(child) != 0 else name)

    return output


def find_missing_chunks(file_name_tree, depth, threshold, file_prefix=""):
    if depth == 0:
        if get_leaf_count(file_name_tree) < threshold:
            print("\t" + file_prefix)
        return
    elif file_name_tree is None or len(file_name_tree) == 0:
        return
    else:
        names = list(file_name_tree.keys())
        names.sort()
        for name in names:
            if name in EXCLUDED_PARTS: continue
            find_missing_chunks(file_name_tree[name], depth-1, threshold,
                                name if file_prefix == "" else file_prefix+"."+name)


def main():
    args = parse_args()
    files = os.listdir(args.input_directory)
    file_name_tree = dict()
    for file in files:
        file = file.strip().split(".")
        add_file_to_name_tree(file_name_tree, file)

    if not args.find_chunks:
        print("File prefixes at depth {}:".format(args.depth))
        for name in get_file_names_for_depth(file_name_tree, args.depth):
            print("\t" + name)
    else:
        print("Missing Chunks at depth {} with threshold {}".format(args.depth, args.threshold))
        find_missing_chunks(file_name_tree, args.depth, args.threshold)


if __name__ == "__main__":
    main()