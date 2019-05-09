#!/usr/bin/env python

from __future__ import print_function

import contextlib
import functools
import os
import re
import sys


def mkdir(dir, mode=0o0750):
    """Construct a directory hierarchy using the given permissions."""
    if not os.path.exists(dir):
        os.makedirs(dir, mode)


@contextlib.contextmanager
def working_directory(path):
    """Changes to the given directory, returning to the original working directory when the context block is exited."""
    original_directory = os.getcwd()
    try:
        os.chdir(path)
        yield
    finally:
        os.chdir(original_directory)


def symlink(source_path, dest_path, abs=False):
    """Create a symbolic link from the source_path to the dest_path, which can be a directory."""

    workdir = os.path.isdir(dest_path) and dest_path or os.path.dirname(dest_path)

    with working_directory(workdir):
        src = abs and os.path.abspath(source_path) or os.path.relpath(source_path, dest_path)
        dest = dest_path
        dest_base = os.path.basename(dest)
        if os.path.isdir(dest_path):
            dest = os.path.join(dest_path, os.path.basename(src))
            if os.path.lexists(dest):
                os.unlink(dest)
            os.symlink(src, dest)
        else:
            mkdir(os.path.dirname(dest_path))
            if os.path.lexists(dest):
                os.unlink(dest)
            os.symlink(src, dest)
        return dest, dest_base
