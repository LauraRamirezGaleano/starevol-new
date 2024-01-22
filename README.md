## FRENCH VERSION

Ce fichier README décrit les spécificités liées à cette branche. Pour une description générale du Git repository, veuillez vous reporter au README de la branche main.

### Nouvelle carte des paramètres (version 4.00)

Cette version de STAREVOL utilise une **nouvelle carte de paramètres (version 4.00) sous format namelist**. Cette nouvelle carte des paramètres simplifie l'ajout et la modification de paramètres et possède une documentation intégrée, basée sur le document starevolpar_340.tex. Les noms des paramètres suivants ont été modifiés par rapport à la version 3.40 :
- paramètre fover ajouté (différentié de etaturb)
- grad_crit &rarr; nu_add
- semiconve &rarr; ledouxmlt
- PM &rarr; PM_turbul
- Dv, Om_max &rarr; Dv_prescr, om_sat
- eq_thermique &rarr; thermal_equilibrium
- Ov, Thl, Tac, Mic &rarr; lover, lthal, ltach, lmicro
- hydro &rarr; hydrodynamics 
- dvarma, dvarmi, dnucma &rarr; dlnvma, dlnvmi, dlnenuc
- itmin &rarr; itermin
- u, r, lnf, lnt, l, ang &rarr; tol_u, tol_lnr, tol_lnf, tol_lnT, tol_l, tol_ang
- resolv &rarr; vlconv

**Pour utiliser la nouvelle carte de paramètres, il vous faut télécharger les fichiers suivants :**
- **starevol.par** : nouvelle carte des paramètres.
- **network.par** : carte des réactions nucléaires, autrefois présente dans starevol.par.
- **rmodpar.f** : nouvelle subroutine de lecture de la carte des paramètres.
- **rininet.f** : nouvelle subroutine de lecture de la carte des réactions nucléaires.
- **starevol.f** : subroutine de STAREVOL initialisant le calcul et qui lit maintenant network.par au lieu de starevol.par.
- **starevol_init.f** : subroutine de STAREVOL initialisant le calcul et qui lit maintenant network.par au lieu de starevol.par.

Le nom de l'ancienne variable **itmin doit également être changé en itermin** dans **nwtraf.f** et dans **evolcom.teq**, celui de **eq_thermique en thermal_equilibrium** dans **prvar.f**, **jtranspc_ondesexcit.f** et **evolcom.diff**, et celui de **grad_crit en nu_add** dans **evolcom.transp ainsi que dans les éventuelles subroutines utilisant grad_crit**. Le paramètre **fover a de plus été ajouté** à la carte des paramètres et est un paramètre global du code, stocké dans un COMMON. **Il doit donc être ajouté dans evolcom.conv dans le COMMON /ovsh/. Les subroutines utilisant fover comme paramètre d'entrée doivent être modifiée pour ne plus l'utiliser comme paramètre d'entrée, et leur appel doit être modifié en conséquence.**

**Les fichiers massloss.par et transport_param.par ne sont plus utilisés. Leurs paramètres ont été ajouté à la nouvelle carte des paramètres.**


## ENGLISH VERSION

This README file describe the specificities of this branch. For a general description of this Git repository, please check the README file of the main branch.

### New parameter card (version 4.00)

This version of STAREVOL uses a **new parameter card (version 4.00) under the namelist format**. This new parameter card simplifies adding and modifing parameters and have an integrated documentation, based on the starevolpar_340.tex document. The names of the following parameters were modified compared to the 3.40 version :
- parameter fover added (differentiated from etaturb)
- grad_crit &rarr; nu_add
- semiconve &rarr; ledouxmlt
- PM &rarr; PM_turbul
- Dv, Om_max &rarr; Dv_prescr, om_sat
- eq_thermique &rarr; thermal_equilibrium
- Ov, Thl, Tac, Mic &rarr; lover, lthal, ltach, lmicro
- hydro &rarr; hydrodynamics 
- dvarma, dvarmi, dnucma &rarr; dlnvma, dlnvmi, dlnenuc
- itmin &rarr; itermin
- u, r, lnf, lnt, l, ang &rarr; tol_u, tol_lnr, tol_lnf, tol_lnT, tol_l, tol_ang
- resolv &rarr; vlconv

**To use the new parameter card, you must download the following files :**
- **starevol.par** : new parameter card.
- **network.par** : nuclear reactions card, before included in starevol.par.
- **rmodpar.f** : new parameter card reading subroutine.
- **rininet.f** : new nuclear reactions card reading subroutine.
- **starevol.f** : STAREVOL subroutine which initialize the computation, which now reads network.par instead of starevol.par.
- **starevol_init.f** : STAREVOL subroutine which initialize the computation, which now reads network.par instead of starevol.par.

The name of the former **itmin variable must also be changed to itermin** in **nwtraf.f** and **evolcom.teq**, that of **eq_thermique to thermal_equilibrium** in **prvar.f**, **jtranspc_ondesexcit.f** and **evolcom.diff**, and that of **grad_crit to nu_add** in **evolcom.transp and in the eventual subroutines using it**. Furthermore **the fover parameter was added** to the parameter card and is now a global parameter of the code, stored in a COMMON. **It must thus be added in evolcom.conv in the /ovsh/ COMMON. The subroutines having fover as an entry parameter must be modified to not have it as an entry parameter anymore, and their call must be modified accordingly.**

**The massloss.par and transport_param.par files are not used anymore. Their parameters were added to the new parameter card.**
