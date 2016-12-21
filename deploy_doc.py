#!/usr/bin/env python
#  -*- coding: utf-8 -*-

from git import Repo
from subprocess import call
import os

from distutils.dir_util import copy_tree, remove_tree
from tempfile import mkdtemp
import sys


def checkout_branch(repo, branch):
    # TODO: faire le check que la branche demandee existe bien dans le repo
    print "\nChecking out to branch %s" % branch
    try:
        print repo.git.checkout(branch)
    except:
        raise RuntimeError('Unable to checkout to branch %s' % branch)


def stash_modifs(repo):
    print "\nStash local modifications"
    print "-------------------------\n"
    try:
        print repo.git.stash(all=True)
    except:
        raise RuntimeError('Unable to stash modifications')


def unstash_modifs(repo):
    print "\nUnstash local modifications"
    print "---------------------------\n"
    try:
        print repo.git.stash('pop', 'stash@{0}')
    except:
        raise RuntimeError('Unable to unstash modifications')


def sphinx_build(working_dir, doc_dir, build_dir='.build'):
    try:
        os.chdir(os.path.join(working_dir, doc_dir))
    except:
        raise RuntimeError('Unable to cd to %s' % doc_dir)
    
    try:
        call(['make', 'clean'])
        call(['make', 'html'])
    except:
        raise RuntimeError('Unable to build Shpinx documentation')
    
    temp_dir = mkdtemp()
    copy_dir_content(os.path.join(working_dir, doc_dir, build_dir, 'html'), temp_dir)
    
    call(['make', 'clean'])
    os.chdir(working_dir)
    
    return temp_dir


def empty_dir():
    pass


def copy_dir_content(src, dst):
    copy_tree(src, dst)
    print "HTML files copied to %s" % dst
    


def retrieve_copy_temp():
    pass


def commit(repo, msg):
    pass


def push_github(repo, remote):
    pass


def get_current_branch(repo):
    branches = repo.git.branch()
    for branch in branches.split('\n'):
        if branch.startswith('*'):
            cur_branch = branch[2:]
            break
    return cur_branch


if __name__ == '__main__':
    
    print("\n========================================")
    print("Deploying Sphinx documentation on GitHub")
    print("========================================")
    
    repo = Repo()
    working_dir = repo.working_tree_dir
    
    # Getting the current branch
    cur_branch = get_current_branch(repo)
    
    stash_modifs(repo)
    
    build_branch = 'master'
    
    try:
        
        # Checking out the master branch
        # TODO: devrait pourvoir etre mis en parametre de la ligne de commande
        checkout_branch(repo, build_branch)
        
        # Building sphinx documentation
        doc_dir = os.path.join(working_dir, 'doc')
        temp_dir = sphinx_build(working_dir, doc_dir)
        
        # Checking out to branch gh-pages
        checkout_branch(repo, 'gh-pages')
        
        # Getting files changed at the last commit
        # last hash
        last_sha = repo.git.log(n='1', pretty="format:%H")
        print last_sha
        # file_list = repo.git.diff_tree(no_commit_id=True, name_only=True, r=last_sha)
        # print file_list
        files = repo.git.diff_tree(no_commit_id=True, name_only=True, r=last_sha).split('\n')
        
        # Removing files
        for file in files:
            try:
                print 'Deleting file %s' % file
                os.remove(file)
            except OSError:
                pass
        
        # Copying new files
        copy_dir_content(temp_dir, '.')
        
        repo.git.add('.', all=True)
        # Commiting changes
        repo.git.commit(m="Updating documentation with respect to branch %s" % build_branch)
        
        
    
    except:
        print "Failed, getting back to initial state"
        checkout_branch(repo, cur_branch)
        unstash_modifs(repo)
        raise RuntimeError('Failed')
    
    checkout_branch(repo, cur_branch)
    unstash_modifs(repo)





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
