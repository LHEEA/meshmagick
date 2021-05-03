
Etude hydrostatique avec meshmagick

On veut ici établir un workflow semi-auto de manière à faire les études hydrostatiques rapidement



* On a un nouveau fichier de maillage

* On vérifie que le maillage est conforme, qu'il n'as pas de pb (-hm), que les normales sont bien orientées

* On vérifie que les dimensions sont bonnes ! Souvent en export Rhino, on se retrouve avec des dimensions en mm, il faut faire un scale -s 0.001

* On vérifie qu'on a bien Z vers le haut

* On vérifie qu'on a un repere qui se trouve:
    - sur le plan de symétrie de la carène OXZ
    - sur la baseline
    - en AP (souvent la mèche de safran)

* Si pas le cas, on doit proposer d"ffectuer les décalages qui vont bien

* On vérifie qu'on a une bonne densité de l'eau

* On s'assure que les données fournies sont bien exprimée:
    - COG dans le repère d'origine du maillage





