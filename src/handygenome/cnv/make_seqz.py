import os
import sys
import re
import subprocess
import multiprocessing
import tempfile

import sequenza.commands

import importlib
top_package_name = __name__.split('.')[0]
common = importlib.import_module('.'.join([top_package_name, 'common']))


def 