## FRENCH VERSION

Voici un tutoriel Git pour pouvoir commencer à collaborer sur ce dépôt. Ce tutoriel est basé sur [celui proposé par Gitlab](https://docs.gitlab.com/ee/gitlab-basics/start-using-git.html).
Prérequis : avoir installé Git sur son ordinateur.

# Cloner le dépôt

Pour cloner ce dépôt, cliquer sur *Clone* et copier l'URL. Saisir 'git clone' dans la ligne de commande suivi de l'URL du dépôt.

# Télécharger les derniers changement apportés au projet

Saisir 'git pull origin <nom-de-la-branche>' en remplaçant par le nom de la branche dans laquelle vous vous trouvez (la branche main est la branche par défaut).

# Changer de branche

Pour aller sur une autre branche, saisir 'git checkout <nom-de-la-branche>. Par exemple, 'git checkout gyre' permet d'accéder à la version du code qui intègre GYRE dans STAREVOL.

# Créer sa propre branche de fonctionnalités

Nous avons opté pour un workflow par branche de fonctionnalités, ce qui signifie que le main contient la version définitive du code et que les nouvelles fonctionnalités sont développées dans des branches dédiées. Pour créer sa propre branche, saisir 'git checkout -b <nom-de-la-branche>'.

# Inspecter les derniers changements

'git diff' permet d'obtenir les différences entre la version locale du projet que vous venez de modifier et la dernière version que vous avez téléchargée via 'git pull' ou 'git clone'.

'git status' répertorie les fichiers qui ont été modifiés, supprimés ou ajoutés.

# Publier des modifications

Plusieurs étapes sont nécessaires :




## ENGLISH VERSION
