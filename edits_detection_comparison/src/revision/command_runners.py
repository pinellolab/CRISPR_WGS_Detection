"""Function running the input commands
"""

import multiprocessing
import subprocess

def run_command(command):
    """Run the input shell command

    :param command: command
    :type command: str
    :return: success
    :rtype: int
    """
    return subprocess.call(command, shell=True)

def run_commands(commands, threads):
    """Run input commands

    :param commands: input commands
    :type commands: List[str]
    :param threads: threads
    :type threads: int
    :raises OSError: raise on command failure
    """
    pool = multiprocessing.Pool(processes=threads)
    result = pool.map_async(run_command, commands)
    codes = result.get()
    for code, cmd in zip(codes, commands):
        if code != 0:
            raise OSError("Edits detction failed!\n\tcommand: %s" % (cmd))
