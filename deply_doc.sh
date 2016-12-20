#!/usr/bin/env bash

cd doc
make clean
make html
cd ..

git add -A
git commit -m "Building docs"

git checkout gh-pages
rm -rf .
touch .nojekyll
git checkout master doc/.build/html
mv ./doc/.build/html/* .
rm -rf ./doc
git add -A
git commit -m "Publishing updated docs"

git checkout master
