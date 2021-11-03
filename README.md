
## FRENCH VERSION

Ce dépôt Git est destiné au développement du code de STAREVOL et à la gestion de ses versions.

- La branche "main" contient la version du code utilisée actuellement. Elle reçoit les changements définitifs du code qui vont être utilisés par toute la communauté.
- La branche "runevol" contient la version de STAREVOL dans laquelle les tables ne sont lues qu'une fois, grâce au programme runevol.f.
- La branche "gyre" est une sous-branche de la branche runevol (donc elle intègre les mêmes changements que la branche "runevol"). Elle contient également le code GYRE adapté pour STAREVOL et le compilateur ifort (dans le dossier gyre-stel), les fichiers d'entrée gyre.in et criteria.in, et les routines modifiées permettant d'appeler GYRE dans le code d'évolution stellaire.
- La branche "freeEOS" est une sous-branche de la branche "gyre". Elle contient les routines permettant d'appeler la librairie FreeEOS.

## ENGLISH VERSION

This Git repository is made for the development and versioning of the STAREVOL code. 

- The main branch contains the version of the code which is currently used. It will receive the final changes of the code which will be used by the whole community.
- The runevol branch contains the STAREVOL version in which tables are only read once, thanks to the runevol.f program.
- The gyre branch is a subbranch of the runevol branch (so it contains the same changes as the runevol branch). It also contains the GYRE code adapted for STAREVOL and the ifort compiler (in the gyre-stel directory), the gyre.in and criteria.in input files, and modified routines allowing to call GYRE in the stellar evolution code.
- The freeEOS branch is a subbranch of the gyre branch. It contains routines calling the FreeEOS library. 

