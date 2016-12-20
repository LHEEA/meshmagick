Ce qui doit être fait pour faire une release open source de Meshmagick
======================================================================


* Migrer sous GitHub:
    * Recupérer des droits sous le compte LHEEA auprès de Guillaume Ducrozet

* Faire l'intégration continue Travis

* Réécrire propremement le README.rst

* Finir la documentation:

    * La partie sur l'hydrostatique

    * Les docstrings des modules

    * Correction des remarques de la doc imprimée


* Uploader la doc en pages perso git en suivant ce `post <http://lucasbardella
  .com/blog/2010/02/hosting-your-sphinx-docs-in-github>`_

* Voir pour un déploiement auto de la doc

* Cleaner les branches du dépôt, faire des rebase voir

* Passer sous `git flow <https://danielkummer.github.io/git-flow-cheatsheet/>`_

* Voir où mettre les conda recipes --> conda forge ?

* Fixer le pb avec programoutput

* Faire que conf.py cherche le numéro de version automatiquement

* Voir d'ailleurs le meilleur endroit pour placer le numéro de version, de manière unique

* S'assurer qu'on a bien la mention GPL partout et qu'on est plus sous la CeCill

* Faire que d-ice beast soit la machine qui fait tourner le code (ou alors la mienne...)

* Ajouter un projet pytriangle et l'option pour créer le lid sur les maillages

* Passer pytriangle sous Git

* Packager pytriangle dans le canal frongere -> packaging avec cython et compilation gcc...

* Créer un environnement dans lequel on est certain de faire tourner meshmagick et le packager dans anaconda
