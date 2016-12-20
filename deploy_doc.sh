#!/usr/bin/env bash

cd doc
make clean
make html
cd ..

git add -A
git commit -m "Building and pushing docs"
git push origin master

git checkout gh-pages
rm -rf .
touch .nojekyll
git checkout master doc/.build/html
mv ./doc/.build/html/* ./
rm -rf ./doc
git add -A
git commit -m "Publishing updated docs"
git push origin gh-pages

git checkout master
