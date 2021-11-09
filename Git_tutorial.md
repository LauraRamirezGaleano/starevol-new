## FRENCH VERSION

Voici un tutoriel Git pour pouvoir commencer à collaborer sur ce dépôt. Ce tutoriel est basé sur [celui proposé par Gitlab](https://docs.gitlab.com/ee/gitlab-basics/start-using-git.html).
Prérequis : avoir installé Git sur son ordinateur.

### Cloner le dépôt

Pour cloner ce dépôt, cliquer sur **Clone** et copier l'URL indiqué dans **Clone with HTTPS**. Saisir `git clone` dans la ligne de commande suivi de l'URL du dépôt.

### Télécharger les derniers changements apportés au projet

Saisir `git pull origin <nom-de-la-branche>` en remplaçant par le nom de la branche dans laquelle on se trouve (la branche main est la branche par défaut).

### Lister les branches 

Plusieurs actions sont possibles :
- savoir sur quelle branche on se situe parmi les branches locales : `git branch` (la branche sur laquelle on se situe est marquée d'un astérisque (*)).
- lister les branches du serveur distant : `git branch -r`.
- lister toutes les branches, locales et distantes : `git branch -a`.

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

### Créer sa propre branche de fonctionnalité

Nous avons opté pour un workflow de branche de fonctionnalité, ce qui signifie que le main contient la version définitive du code et que les nouvelles fonctionnalités sont développées dans des branches dédiées. Pour créer sa propre branche, il est possible de cliquer sur **+** puis **New branch** depuis la page GitLab. Sinon, pour créér la branche depuis la ligne de commande :
- se placer sur le main : `git checkout main`. 
- télécharger les derniers changements apportés au main : `git pull origin main`.
- créer sa branche : `git checkout -b <nom-de-la-nouvelle-branche>`.

Une fois la branche créée, ne pas oublier de publier les changements sur le dépôt distant : `git push origin <nom-de-la-nouvelle-branche>`.

### Fusionner une branche avec le main

Cette action peut être effectuée lorsque l'on veut intégrer une fonctionnalité développée sur une branche à la version définitive du code qui se trouve sur le main. Pour cela : 
- se placer sur le main : `git checkout main`.
- télécharger les derniers changements apportés au main : `git pull origin main`. 
- fusionner la branche avec le main : `git merge <nom-de-la-branche>`.

### Intégrer des changements d'autres branches dans sa branche

Pour intégrer des changements d'autres branches (par exemple les changements de la branche gyre) :
- télécharger les derniers changements apportés à la branche gyre : `git checkout gyre` puis `git pull origin gyre`.
- se placer sur sa branche : `git checkout <nom-de-sa-branche>`.
- fusionner la branche gyre avec sa branche : `git merge gyre`.

## ENGLISH VERSION

This it a Git tutorial to start collaborating on this repository. This tutorial is based on [the GitLab one](https://docs.gitlab.com/ee/gitlab-basics/start-using-git.html).
Prerequisite : having Git installed on your computer.

### Clone the repository

To clone the repository, click on **Clone** and copy the URL written in **Clone with HTTPS**. Run `git clone` in the command line and paste the repository URL.

### Download the latest changes

Run `git pull origin <branch-name>` (replace by the branch name in which you are - main is the default branch).

### List branches 

Several actions are possible :
- to know on which branch you are among local branches : `git branch` (the current local branch will be marked with an asterisk (*)).
- to list remote branches : `git branch -r`.
- to list all remote and local branches : `git branch -a`.

### Switch to a branch

To switch to another branch, run `git checkout <branch-name>`. For example, `git checkout gyre` allows you to go to the code version which contains GYRE for STAREVOL.

### View latest changes

`git diff` allows you to view the differences between your local version of the project which you are modifying and the latest version you downloaded via `git pull` or `git clone`.

`git status` allows you to view which files have been modified, removed or added.

### Push modifications

You will need to follow several steps :

- add changes : `git add <file-name or directory-name>`. For example, if you are in the git directory : `git add .` will add every changed made in this directory.
- commit changes : `git commit -m "Commit message"` or `git commit` if you want to use a text editor to write the commit message (which is convenient for long messages).
- push changes on the remote git repository : `git push origin <branch-name>`.

### Create your own feature branch

We agreed on a feature branch workflow, which means that the main branch contains the final version of the code and the new features are developed in dedicated branches. To create your own branch, you can click on **+** then **New branch** on the GitLab page. Otherwise, you can create the branch from the command line :
- switch to the main : `git checkout main`.
- download the latest changes of the main : `git pull origin main`. 
- create your branch : `git checkout -b <new-branch-name>`.

Once your branch is created, do not forget to push your changes to the remote repository : `git push origin <new-branch-name>`.

### Merge a branch with the main

If you want to integrate a feature developed in a branch to the final version of the code which is on the main, you need to :
- switch to the main : `git checkout main`.
- download the latest changes of the main : `git pull origin main`. 
- merge the feature branch to the main : `git merge <branch-name>`.

### Get features from other branches in your branch

To get features from other branches (for example the gyre branch features) :
- download the latest changes of the gyre branch : `git checkout gyre` then `git pull origin gyre`. 
- switch to your own branch : `git checkout <your-branch-name>`.
- merge the gyre branch to your own branch : `git merge gyre`.


