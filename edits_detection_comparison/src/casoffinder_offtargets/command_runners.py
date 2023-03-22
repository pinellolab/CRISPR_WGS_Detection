from tqdm import tqdm

import multiprocessing
import subprocess

def run_command(command):
    """Run the command

    :param command: command
    :type command: str
    :return: command signal
    :rtype: int
    """
    return subprocess.call(command, shell=True)


def run_commands(commands, threads):
    """Run the input commands

    :param commands: input commands
    :type commands: List[str]
    :param threads: threads
    :type threads: int
    :raises OSError: raise on Mutect2 failure
    """
    pool = multiprocessing.Pool(processes=threads)  # create `threads` threads
    result = pool.map_async(run_command, commands)  # run commands in parallel
    with tqdm(total=len(commands)) as progress:
        while not result.ready():
            remaining = result._number_left
            progress.update(len(commands) - remaining)
    codes = result.get()
    for code, cmd in zip(codes, commands):
        if code != 0:
            raise OSError("An error occurred while running %s" % (cmd))