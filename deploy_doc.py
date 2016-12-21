#!/usr/bin/env python
#  -*- coding: utf-8 -*-

from git import Repo
import os
from os import chdir, system, getcwd
from distutils.dir_util import copy_tree, remove_tree
from tempfile import mkdtemp
import sys


def checkout_branch(repo, branch):
    try:
        repo.git.checkout(branch)
    except:
        raise RuntimeError('Unable to checkout to branch %s' % branch)


def stash_modifs(repo):
    try:
        repo.git.stash()
    except:
        raise RuntimeError('Unable to stash modifications')
    
    
def unstash_modifs(repo):
    try:
        repo.git.stash('pop')
    except:
        raise RuntimeError('Unable to unstash modifications')


def sphinx_build(repo):
    pass


def empty_dir():
    pass


def copy_temp():
    return 'path'


def retrieve_copy_temp():
    pass


def commit(repo, msg):
    pass


def push_github(repo, remote):
    pass


if __name__ == '__main__':
    
    print("Deploying documentation")
    print("-----------------------\n")
    
    repo = Repo()
    
    # Stash modifications
    stash_modifs(repo)
    
    
    



    
# repo = Repo()
# git = repo.git
#
# branches = git.branch()
# print branches
#
# branches.split()
#
# # Getting the current branch (better way to do that ?)
# for branch in branches.split('\n'):
#     if branch.startswith('*'):
#         cur_branch = branch[2:]
#
#
# # Checking out to master branch
# git.checkout('master')
#
# # Going to doc dir and building documentation
# chdir('doc')
# system('make clean')
# system('make html')
#
# # Copying html files into a tempdir
# tempdir = mkdtemp()
# print tempdir
#
# # TODO: aller cherche dans le conf.py le dossier build (.build ou _build)
# html_files = os.path.join(getcwd(), '.build/html')
#
#
# copy_tree(html_files, tempdir)
# chdir('..')
#
#
# git.checkout('gh-pages')
#
# # Removing everything here
# system('rm -rf ./*')  # TODO: remplacer par des fonctions pures python
#
# copy_tree(tempdir, '.')
#
# git.commit(m='Updating documentation', a=True)
#
#
# # Checking out back to the current branch
# git.checkout(cur_branch)
#
#
# print git.branch()
