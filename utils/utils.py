import os
import subprocess


def rel_to_abs(rel_path, file):
    """
    file is __file__ of wherever relpath starts
    """
    return os.path.join(os.path.dirname(file), rel_path)


def system_call(command):
    p = subprocess.Popen([command], stdout=subprocess.PIPE, shell=True)
    return p.stdout.read()