## FRENCH VERSION

Voici un tutoriel Git pour pouvoir commencer à collaborer sur ce dépôt. Ce tutoriel est basé sur [celui proposé par Gitlab](https://docs.gitlab.com/ee/gitlab-basics/start-using-git.html).
Prérequis : avoir installé Git sur son ordinateur.

### Cloner le dépôt

Pour cloner ce dépôt, cliquer sur **Clone** et copier l'URL. Saisir `git clone` dans la ligne de commande suivi de l'URL du dépôt.

### Télécharger les derniers changement apportés au projet

Saisir `git pull origin <nom-de-la-branche>` en remplaçant par le nom de la branche dans laquelle on se trouve (la branche main est la branche par défaut).

### Changer de branche

Pour aller sur une autre branche, saisir `git checkout <nom-de-la-branche>`. Par exemple, `git checkout gyre` permet d'accéder à la version du code qui intègre GYRE dans STAREVOL.

### Inspecter les derniers changements

`git diff` permet d'obtenir les différences entre la version locale du projet que l'on vient de modifier et la dernière version téléchargée via `git pull` ou `git clone`.

`git status` répertorie les fichiers qui ont été modifiés, supprimés ou ajoutés.

### Publier des modifications

Plusieurs étapes sont nécessaires :

- indexer les changements : `git add <nom-du-fichier ou du dossier>`. Par exemple, si l'on se trouve dans le dossier git : `git add .` permet d'indexer tous les changements effectués dans ce dossier.

- enregistrer les changements : `git commit -m "Message de commit"` ou `git commit` si l'on veut taper le message de commit dans un éditeur de texte (pratique pour les messages détaillés).

- publier les changements sur le dépôt du projet : `git push origin <nom-de-la-branche>`.

### Fusionner une branche avec le main

- Se placer sur le main : `git checkout main`.
- Fusionner la branche avec le main : `git merge <nom-de-la-branche>`.


### Créer sa propre branche de fonctionnalité

Nous avons opté pour un workflow de branche de fonctionnalité, ce qui signifie que le main contient la version définitive du code et que les nouvelles fonctionnalités sont développées dans des branches dédiées. Pour créer sa propre branche :
- se placer sur le main : `git checkout main`. 
- créer sa branche : `git checkout -b <nom-de-la-nouvelle-branche>`.

Pour intégrer des changements d'autres branches (par exemple les changements de la branche gyre) :
- se placer sur sa branche : `git checkout <nom-de-la-nouvelle-branche>`.
- fusionner la branche gyre avec sa branche : `git merge gyre`.


## ENGLISH VERSION

This it a Git tutorial to start collaborating on this repository. This tutorial is based on [the GitLab one](https://docs.gitlab.com/ee/gitlab-basics/start-using-git.html).
Prerequisite : having Git installed on your computer.

### Clone the repository

To clone the repository, click on **Clone** and copy URL. Run `git clone` in the command line and paste the repository URL.

### Download the latest changes

Run `git pull origin <branch-name>` (replace by the branch name in which you are - main is the default branch).

### Switch to a branch

To switch to another branche, run `git checkout <branch-name>`. For example, `git checkout gyre` allows you to go the code version which contains GYRE for STAREVOL.

### View latest changes

`git diff` allows you to view the differences between your local version of the project which you are modifying and the latest version you downloaded via `git pull` or `git clone`.

`git status` allows you to view which files have been modified, removed or added.

### Push modifications

You will need to follow several steps :

- add changes : `git add <file-name or directory-name>`. For example, if you are in the git directory : `git add .` will add every changed made in this directory.

- commit changes : `git commit -m "Commit message"` or `git commit` if you want to use a text editor to write the commit message (which is convenient for long messages).

- push changes on the remote git repository : `git push origin <branch-name>`.

### Merge a branch with the main

- Switch to the main : `git checkout main`.
- Merge the feature branch to the main : `git merge <branch-name>`.


### Create your own feature branch

We agreed on a feature branch workflow, which means that the main branch contains the final version of the code and the new features are developed in dedicated branches. To create your own branch :
- switch to the main : `git checkout main`.
- create your branch : `git checkout -b <new-branch-name>`.

To get features from other branches (for example the gyre branch features) :
- switch to your own branch : `git checkout <new-branch-name>`.
- merge the gyre branch to your own branch : `git merge gyre`.
