#!/usr/bin/env python
#  -*- coding: utf-8 -*-

from git import Repo
import os
from os import chdir, system, getcwd
from distutils.dir_util import copy_tree, remove_tree
from tempfile import mkdtemp
import sys

repo = Repo()
git = repo.git

branches = git.branch()
print branches

branches.split()

# Getting the current branch (better way to do that ?)
for branch in branches.split('\n'):
    if branch.startswith('*'):
        cur_branch = branch[2:]


# Retier
git.commit(m='Update deploy', a=True)


sys.exit(0)

# Checking out to master branch
git.checkout('master')

# Going to doc dir and building documentation
chdir('doc')
system('make clean')
system('make html')

# Copying html files into a tempdir
tempdir = mkdtemp()
print tempdir

# TODO: aller cherche dans le conf.py le dossier build (.build ou _build)
html_files = os.path.join(getcwd(), '.build/html')


copy_tree(html_files, tempdir)
chdir('..')


git.checkout('gh-pages')

# Removing everything here
print getcwd()


system('rm -rf ./*')  # TODO: remplacer par des fonctions pures python

copy_tree(tempdir, '.')


git.commit("")
# Commit de tout ca


sys.exit(0)





# Checking out back to the current branch
git.checkout(cur_branch)


print git.branch()
