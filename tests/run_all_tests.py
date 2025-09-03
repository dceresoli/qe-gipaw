#!/usr/bin/env python
# run all GIPAW tests at once, comparing results with the previous version

import sys, os, subprocess
import shlex

# helper function to run a shell command
def run(cmd, stdin=None, stdout=None, stderr=None, echo=False):
    if echo:
        echo_msg = 'running: ' + cmd
        if stdin is not None:
            echo_msg += ' <' + stdin
        if stdout is not None:
            echo_msg += ' >' + stdout
        if stderr is not None:
            echo_msg += ' 2>' + stderr
        print(echo_msg)

    if stdin is not None:
        stdin = open(stdin, 'rt')
    if stdout is not None:
        stdout = open(stdout, 'wt')

    with subprocess.Popen(shlex.split(cmd), bufsize=1024, stdin=stdin, stdout=stdout) as proc:
        proc.wait()

    if stdin is not None:
        stdin.close()
    if stdout is not None:
        stdout.close()

    return proc.returncode


# EFG
def test_EFG():
    global VERSION, QEDIR, testdir

    for subdir in ['NC', 'US', 'PAW']:
        print(f'testing EFG/{subdir}...')
        os.chdir(f'{testdir}/EFG/{subdir}')
        run('bash run_tests.sh')
    print()


# NMR
def test_NMR():
    global VERSION, QEDIR, testdir

    for subdir in ['NC', 'US', 'PAW']:
        print(f'testing NMR/{subdir}...')
        os.chdir(f'{testdir}/NMR/{subdir}')
        run('bash run_tests.sh')
    print()


# Hubbard
def test_Hubbard():
    global VERSION, QEDIR, testdir

    print(f'testing Hubbard...')
    os.chdir(f'{testdir}/Hubbard')
    run('bash run_tests.sh')
    print()


# EPR
def test_EPR():
    global VERSION, QEDIR, testdir

    print(f'testing EPR...')
    os.chdir(f'{testdir}/EPR')
    run('bash run_tests.sh')
    print()



# get the current version
with open('version', 'r') as f:
    exec(f.read())

print('='*72)
print(f'QE version: {VERSION}, found in {QEDIR}')
print('='*72)
print()

# save the current directory 
testdir = os.getcwd()

test_EFG()
test_NMR()
test_Hubbard()
test_EPR()
