# Analyse des Variants Génétiques de la Cardiomyopathie Hypertrophique dans ClinVar

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![Bash 4.4+](https://img.shields.io/badge/bash-4.4+-green.svg)](https://www.gnu.org/software/bash/)
[![Dataset ClinVar](https://img.shields.io/badge/data-NCBI%20ClinVar-orange.svg)](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/archive/variant_summary_2026-10.txt.gz)
[![License: CC0 / Public Domain](https://img.shields.io/badge/license-Public%20Domain-lightgrey.svg)](SOURCES.md)

Ce projet implémente un pipeline bioinformatique automatisé et reproductible combinant **Bash** et **Python 3** pour identifier et hiérarchiser les gènes concentrant le plus grand nombre de variants distincts classés **« Pathogenic »** ou **« Likely pathogenic »** dans la base de données **NCBI ClinVar**, spécifiquement associés à la **cardiomyopathie hypertrophique (CMH / HCM)**.

> **Groupe 11** · Dépôt d'équipe : <https://github.com/anfalda/projectSE> (le hash du commit final est dans `LIEN_GIT.txt`).
> Ce README suffit pour installer, lancer et vérifier le pipeline sur une Ubuntu vierge, sans intervention manuelle en dehors des commandes ci-dessous.

---

## Sommaire

1. [Problématique Scientifique et Objet d'Étude](#1-problématique-scientifique-et-objet-détude)
2. [Architecture du Projet](#2-architecture-du-projet)
3. [Chaîne de Traitement (Pipeline)](#3-chaîne-de-traitement-pipeline)
4. [Prérequis et Installation](#4-prérequis-et-installation)
5. [Guide d'Exécution](#5-guide-dexécution)
   - [Exécution complète avec run.sh](#exécution-complète-avec-runsh)
   - [Options de ligne de commande (getopts)](#options-de-ligne-de-commande-getopts)
   - [Exécution modulaire (Filtre Shell & Analyse Python)](#exécution-modulaire-filtre-shell--analyse-python)
6. [Tests Automatiques](#6-tests-automatiques)
7. [Fichiers de Résultats et Interprétation Biologique](#7-fichiers-de-résultats-et-interprétation-biologique)
   - [Classement des principaux gènes](#classement-des-principaux-gènes)
   - [Description des sorties](#description-des-sorties)
8. [Administration Système et Reproductibilité](#8-administration-système-et-reproductibilité)
9. [Dépannage et FAQ](#9-dépannage-et-faq)
10. [Corrections issues du test croisé](#10-corrections-issues-du-test-croisé)
11. [Membres de l'équipe et contributions](#membres-de-léquipe-et-contributions)

---

## 1. Problématique Scientifique et Objet d'Étude

- **Question de recherche :** *Quels gènes concentrent le plus de variants distincts classés « Pathogenic » ou « Likely pathogenic » dans ClinVar, associés spécifiquement à la cardiomyopathie hypertrophique ?*
- **Pathologie étudiée :** La **cardiomyopathie hypertrophique (CMH)** est une maladie myocardique génétique fréquente (prévalence estimée à 1/500 dans la population générale), caractérisée par une hypertrophie ventriculaire gauche inexpliquée par une surcharge de pression ou de volume. Elle constitue l'une des causes majeures de mort subite chez l'adulte jeune et l'athlète.
- **Distinction scientifique essentielle :** Pour chaque gène, le classement compte les `VariationID` distincts. Le total global correspond aux **associations variant–gène retenues** après dédoublonnage au sein de chaque gène. Un même `VariationID` peut être associé à plusieurs gènes. Les pourcentages représentent la part de chaque gène dans le total des associations retenues, et non une prévalence clinique ou un nombre de patients porteurs.

---

## 2. Architecture du Projet

```text
projectSE/
├── assets/
│   └── usthb_logo.png
├── config/
│   └── dataset.env
├── cron/
│   └── clinvar_pipeline.cron
├── data/
│   ├── raw/                     # Données brutes téléchargées par run.sh (non versionnées, lecture seule)
│   └── processed/               # Données filtrées
├── logs/                        # Journaux horodatés (run_*.log)
├── results/                     # Résultats CSV, JSON, PNG et archive datée (ceux de la dernière exécution sont versionnés)
├── scripts/
│   ├── analyze_genes.py         # Analyse Python (argparse, CSV/JSON, figure)
│   └── filter_hcm_variants.sh   # Filtrage Shell (tubes awk | sort | uniq)
├── tests/
│   ├── test_analyze_genes.py
│   └── test_results_documentation.py
├── tmp/                         # Fichiers temporaires (nettoyés par trap)
├── .venv/                       # Environnement virtuel (créé par run.sh, non versionné)
├── .gitattributes
├── .gitignore
├── LIEN_GIT.txt                 # Adresse du dépôt et commit final
├── README.md
├── requirements.txt             # Dépendances Python (matplotlib)
├── rapport.pdf
├── run.sh                       # Point d'entrée unique du pipeline
└── SOURCES.md                   # Identification et intégrité du jeu de données
```

Les données ne sont **pas** fournies avec le projet : `run.sh` télécharge lui-même le fichier ClinVar et vérifie son intégrité.

---

## 3. Chaîne de Traitement (Pipeline)

Le projet applique une séparation claire des responsabilités entre **Bash** (flux système, gestion de fichiers, filtrage de gros volumes) et **Python** (normalisation, structures relationnelles, statistiques et visualisation graphique) :

```
[NCBI ClinVar] -> variant_summary_2026-10.txt.gz (archive 2026-10, ≈ 450,7 MB compressés, soit ≈ 430 MiB)
       │
       ▼ (1. Téléchargement & Contrôle d'intégrité)
  curl / sha256sum / gzip -t
       │
       ▼ (2. Filtrage Shell haute performance en flux)
  scripts/filter_hcm_variants.sh
  - Assemblage : GRCh38
  - Origine : germline
  - Phénotype : Hypertrophic cardiomyopathy (strict)
       │
       ▼
  data/processed/hcm_variants_grch38.tsv (17 012 enregistrements)
       │
       ▼ (3. Analyse Python, dédoublonnage et tri)
  scripts/analyze_genes.py
  - Retient : Pathogenic, Likely pathogenic, Pathogenic/Likely pathogenic
  - Exclut : VUS, Benign, Conflicting, gènes manquants ou multi-gènes (;)
  - Dédoublonne les VariationID au sein de chaque gène
       │
       ▼ (4. Génération des résultats)
  results/gene_ranking.csv, variants_hcm.csv, excluded_variants.csv, summary.json, top_genes.png
       │
       ▼ (5. Archivage horodaté & Nettoyage)
  results/resultats_YYYYMMDD_HHMMSS.tar.gz & trap cleanup
```

**Enchaînement par tubes (étape 2) :** le filtrage lit l'archive en flux, sans jamais décompresser les ~4 Go sur le disque :

```bash
gzip -cd variant_summary_2026-10.txt.gz \
  | awk -F '\t' '...filtre GRCh38 / germline / phénotype HCM...' \
  | sort -t $'\t' -k1,1n -k11,11n \
  | uniq > fichier_filtre.tsv
```

---

## 4. Prérequis et Installation

### Système d'exploitation
- **Linux** (Debian, Ubuntu, RHEL, Arch) ou **WSL 2** (Windows Subsystem for Linux avec Ubuntu).
- Testé sous **Ubuntu 24.04**. Git Bash (Windows) n'est pas garanti, voir la FAQ.

### Outils système requis
- `bash` (version >= 4.4)
- `curl` (téléchargement réseau)
- `gzip` (compression/décompression)
- `awk`, `sort`, `uniq`, `tar`, `tee`, `mktemp`
- `sha256sum` ou `shasum`

### Ressources nécessaires
- Connexion Internet (téléchargement d'environ 430 MiB depuis le FTP du NCBI).
- Environ 1 Go d'espace disque libre (archive brute + environnement virtuel + résultats).
- Durée : environ 55 secondes pour le traitement (hors téléchargement) sur Ubuntu 24.04 (VMware).

### Installation des paquets système (Ubuntu / Debian)
Sur une Ubuntu vierge, `python3-venv` n'est pas installé par défaut : sans lui, la création de l'environnement virtuel échoue.
```bash
sudo apt update
sudo apt install -y git python3 python3-venv python3-pip curl gzip
```

### Environnement Python
- Python 3.10 ou version supérieure.
- Installation dans un environnement virtuel Python :

> **Cron :** le chemin `$HOME/projectSE` est codé en dur dans `cron/clinvar_pipeline.cron`. Le dépôt doit donc être cloné dans un dossier nommé exactement `projectSE`, sinon la tâche planifiée échoue. Si vous utilisez un autre nom, remplacez `projectSE` dans le fichier `.cron` avant de l'installer (voir section 8).

```bash
git clone https://github.com/anfalda/projectSE.git "$HOME/projectSE"
cd "$HOME/projectSE"
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
chmod +x run.sh scripts/*.sh scripts/*.py
./run.sh
```

Effectuer ce premier lancement avant d'activer Cron : il prépare les dossiers et l'environnement.

**Sortie attendue (extrait) :**
```text
[INFO] Demarrage du pipeline ClinVar 2026-10
[OK] Integrite verifiee: variant_summary_2026-10.txt.gz
Lignes HCM strictes conservees: 17012
     Associations retenues    : 1427
     Enregistrements exclus   : 15585
     Genes distincts   : 40
[OK] Pipeline termine avec succes.
```
Un lancement réussi se termine avec le code de retour `0` (`echo $?`).

---

## 5. Guide d'Exécution

### Rendre les scripts exécutables
Avant la première utilisation, assurez-vous que les permissions d'exécution sont positionnées :
```bash
chmod +x run.sh scripts/*.sh scripts/*.py
```

### Exécution complète avec `run.sh`
Pour exécuter l'ensemble de la chaîne de bout en bout (téléchargement conditionnel, filtrage, analyse Python, archivage et logs) :

```bash
./run.sh
```

### Options de ligne de commande (`getopts`)
Le script `run.sh` intègre plusieurs options documentées :

| Option | Argument | Description | Valeur par défaut |
| :--- | :--- | :--- | :--- |
| `-c` | `FILE` | Chemin vers un fichier de configuration alternatif | `config/dataset.env` |
| `-t` | `N` | Nombre de gènes affichés dans la figure finale (entier strictement positif, sinon erreur avec code retour 1) | `15` |
| `-f` | Aucun | Force le re-téléchargement même si le fichier brut existe | Désactivé |
| `-n` | Aucun | Mode *dry-run* : affiche le plan d'exécution sans traitement | Désactivé |
| `-h` | Aucun | Affiche l'aide détaillée et quitte | — |

**Exemples :**
```bash
# Afficher l'aide
./run.sh -h

# Afficher le plan d'exécution sans toucher aux données
./run.sh -n

# Exécuter avec un Top 20 de gènes
./run.sh -t 20

# Forcer le re-téléchargement de la base ClinVar
./run.sh -f
```

### Exécution modulaire (Filtre Shell & Analyse Python)

Si vous souhaitez exécuter séparément chaque composant (l'étape Python nécessite matplotlib, donc l'environnement virtuel doit être activé) :

```bash
source .venv/bin/activate
```

```bash
# 1. Étape Shell seule
./scripts/filter_hcm_variants.sh data/raw/variant_summary_2026-10.txt.gz data/processed/hcm_variants_grch38.tsv

# 2. Étape Python seule
python3 scripts/analyze_genes.py \
  --variants data/processed/hcm_variants_grch38.tsv \
  --output-dir results \
  --top 15
```

---

## 6. Tests Automatiques

Une suite de tests automatisés basée sur le module standard `unittest` de Python est répartie dans `tests/test_analyze_genes.py` et `tests/test_results_documentation.py`. Elle vérifie sur des jeux de données synthétiques indépendants :
1. L'acceptation stricte des classes `Pathogenic`, `Likely pathogenic` et `Pathogenic/Likely pathogenic`.
2. Le rejet des classifications bénignes, incertaines (VUS) et conflictuelles.
3. Le dédoublonnage des `VariationID` identiques pour un même gène.
4. L'exclusion documentée des symboles de gènes manquants (`-`, vide, `NA`).
5. L'exclusion documentée des symboles de gènes ambigus multi-gènes (séparateur `;`).
6. Le tri décroissant exact du classement des gènes.
7. La conformité de l'interface CLI et la création des 5 livrables.
8. La cohérence entre les résultats générés et les statistiques documentées dans le README.

Pour lancer les tests, **activer d'abord l'environnement virtuel** (créé lors du premier `./run.sh` ou à l'étape d'installation). Sans cela, les tests échouent avec `ModuleNotFoundError: No module named 'matplotlib'`.
```bash
cd "$HOME/projectSE"
source .venv/bin/activate
python -m unittest discover -s tests -v
```

---

## 7. Fichiers de Résultats et Interprétation Biologique

### Classement des principaux gènes

Sur le jeu de données ClinVar (Release 2026-10), le filtrage a produit **17 012 enregistrements HCM**. Après les critères de sélection et le dédoublonnage, **1 427 associations variant–gène** sont retenues et réparties sur **40 gènes**. Au total, **15 585 enregistrements** sont exclus.

Les nombres de la colonne « VariationID distincts par gène » sont des décomptes distincts au sein de chaque gène. Les pourcentages sont calculés sur les 1 427 associations variant–gène retenues.

| Rang | Symbole Gène | Nom de la Protéine | VariationID distincts par gène | Part (%) | Rôle Physiopathologique |
| :---: | :--- | :--- | :---: | :---: | :--- |
| **1** | **MYBPC3** | Myosin-binding protein C, cardiaque | **832** | **58,30 %** | Protéine sarcomérique clé ; mutations tronquantes fréquentes (haplo-insuffisance). |
| **2** | **MYH7** | Chaîne lourde de bêta-myosine cardiaque | **309** | **21,65 %** | Protéine du filament épais ; mutations faux-sens dominantes au niveau de la tête motrice. |
| **3** | **FLNC** | Filamine C | **66** | **4,63 %** | Protéine du disque Z et du cytosquelette intercalaire ; associée aux formes arythmogènes. |
| **4** | **TNNI3** | Troponine I cardiaque | **56** | **3,92 %** | Composante du complexe régulateur de la troponine du filament fin. |
| **5** | **TNNT2** | Troponine T2 cardiaque | **32** | **2,24 %** | Régulation calcique de la contraction ; associée à des risques de mort subite précoce. |
| **6** | **TPM1** | Tropomyosine 1 (alpha) | **24** | **1,68 %** | Filament fin du sarcomère. |
| **7** | **TTN** | Titine | **14** | **0,98 %** | Protéine géante de cohésion élastique du sarcomère. |
| **8** | **PRKAG2** | Sous-unité gamma-2 de la protéine kinase AMPK | **11** | **0,77 %** | Phénocopie métabolique de CMH avec préexcitation ventriculaire (WPW). |
| **9** | **MYL2** | Chaîne légère régulatrice de myosine | **10** | **0,70 %** | Filament épais. |
| **10** | **ALPK3** | Alpha-kinase 3 | **8** | **0,56 %** | Régulateur de l'intégrité sarcomérique (formes néonatales/sévères). |

Le tableau ci-dessus présente les 10 premiers gènes ; le classement complet des 40 gènes est dans `results/gene_ranking.csv`, et la figure `top_genes.png` affiche le Top N choisi avec `-t` (15 par défaut).

> **Constat biologique majeur :** MYBPC3 et MYH7 totalisent **79,96 %** des associations variant–gène retenues dans cette analyse. Ce résultat décrit la répartition des enregistrements ClinVar filtrés ; il ne mesure pas la fréquence de ces variants dans la population ni le nombre de patients atteints.

### Description des sorties

- **`results/gene_ranking.csv`** : Tableau récapitulatif classé contenant le rang, le symbole HGNC, le décompte de variants et le pourcentage.
- **`results/variants_hcm.csv`** : Tableau complet des 1 427 associations variant–gène retenues avec leurs identifiants ClinVar (`VariationID`, `AlleleID`, `RCVaccession`, etc.).
- **`results/excluded_variants.csv`** : Traçabilité des 15 585 enregistrements exclus avec le motif précis (classification non retenue, gène ambigu, doublon).
- **`results/summary.json`** : Fichier JSON structuré pour intégration web ou programmatique contenant les métadonnées, statistiques globales et note de portée scientifique.
- **`results/top_genes.png`** : Graphique haute résolution (300 DPI) à barres horizontales prêt pour l'intégration dans le rapport scientifique.
- **`results/resultats_YYYYMMDD_HHMMSS.tar.gz`** : Archive compressée d'intégrité regroupant les 5 fichiers de résultats générés lors de l'exécution.

---

## 8. Administration Système et Reproductibilité

1. **Intégrité cryptographique (SHA-256) :**  
   Avant tout traitement, `run.sh` vérifie l'empreinte SHA-256 du jeu de données brut téléchargé afin de détecter une corruption ou une modification du fichier.
2. **Nettoyage automatique par signaux (`trap`) :**  
   Le script configure `trap cleanup EXIT` et `trap on_signal INT TERM HUP`. Tout fichier temporaire dans `tmp/` est supprimé proprement en cas de succès ou d'interruption par `Ctrl+C`.
3. **Protection des données brutes en lecture seule :**  
   Dès son téléchargement et sa vérification, le fichier `data/raw/variant_summary_2026-10.txt.gz` est verrouillé en lecture seule via `chmod a-w`.
4. **Journaux d'exécution horodatés :**  
   Chaque exécution génère un fichier journal dédié `logs/run_YYYYMMDD_HHMMSS.log` capturant la sortie standard et d'erreur avec un horodatage à la seconde.
5. **Répétabilité sans intervention manuelle :**  
   Le pipeline peut être relancé sans nettoyage manuel préalable ; chaque exécution produit un nouveau journal et une nouvelle archive horodatée.
6. **Planification automatique (Cron) :**  
   Le fichier `cron/clinvar_pipeline.cron` permet de planifier une réexécution mensuelle de l'analyse sur la release archivée `2026-10`. Cette planification assure la répétabilité du traitement, mais ne met pas automatiquement à jour la version du jeu de données.
   - **Chemin codé en dur :** le fichier utilise `$HOME/projectSE`. Si le dépôt porte un autre nom, remplacer `projectSE` dans le fichier avant de l'installer (`sed 's/projectSE/NOM_DU_DOSSIER/g' cron/clinvar_pipeline.cron | crontab -`, en veillant à ne pas écraser un crontab existant).
   - **Nettoyage automatique :** chaque dimanche à 03h00, les journaux `run_*.log` et les archives `resultats_*.tar.gz` de plus de 90 jours sont supprimés.
   - **Prérequis :** lancer `./run.sh` une première fois à la main (création de `.venv` et de `logs/`) avant d'activer la tâche.
   - **Installation de la tâche** (ajoute les lignes sans écraser le crontab existant) :
     ```bash
     (crontab -l 2>/dev/null; cat cron/clinvar_pipeline.cron) | crontab -
     crontab -l        # vérification
     ```
     Planification : réexécution le 1er de chaque mois à 02h00, journal dans `logs/cron_AAAAMMJJ.log`.
7. **Droits d'accès :**
   Les scripts sont exécutables (`chmod +x`), et les données brutes sont en lecture seule (`chmod a-w`). Le script retire temporairement cette protection uniquement pour un re-téléchargement forcé (`-f`).
8. **Vérification manuelle de l'intégrité :**
   L'empreinte attendue est dans `config/dataset.env` (`VARIANT_SHA256`). Pour la recalculer : `sha256sum data/raw/variant_summary_2026-10.txt.gz`.
9. **Codes de retour et messages d'erreur :**

   | Code | Signification |
   | :---: | :--- |
   | `0` | Succès |
   | `1` | Erreur détectée avec message clair (option invalide, `-t` non positif, fichier manquant ou corrompu, commande ou Python introuvable, SHA-256 incorrecte) |
   | `2` | `filter_hcm_variants.sh` : mauvais nombre d'arguments |
   | `3` / `4` | `filter_hcm_variants.sh` : colonne requise manquante / en-tête ClinVar introuvable |
   | `130` | Interruption par signal (`Ctrl+C`, `TERM`, `HUP`) |

   Le mode `set -Eeuo pipefail` arrête le pipeline à la première erreur.
10. **Fonctions réutilisables dans `run.sh` :**
   `log`, `die`, `cleanup`, `on_signal`, `require_command`, `sha256_of`, `verify_file`, `download_file`, `run_shell_filter`, `run_python_analysis`, `archive_results`.

---

## 9. Dépannage et FAQ

- **Erreur `Interpreteur Python introuvable` :**  
  Spécifiez explicitement le binaire Python via la variable `PYTHON_BIN` :  
  `PYTHON_BIN=/usr/bin/python3 ./run.sh`
- **Erreur `Somme SHA-256 incorrecte` :**  
  Vérifiez que le téléchargement n'a pas été tronqué. Utilisez l'option `-f` pour forcer un re-téléchargement complet :  
  `./run.sh -f`
- **Erreur `venv impossible (installez python3-venv)` :**  
  Installer le paquet manquant : `sudo apt install python3-venv python3-pip`, puis relancer `./run.sh`.
- **Erreur `ModuleNotFoundError: No module named 'matplotlib'` :**  
  L'environnement virtuel n'est pas activé : `source .venv/bin/activate` (ou lancer `./run.sh`, qui utilise automatiquement `.venv`).
- **La tâche Cron ne s'exécute pas :**  
  Vérifier que le dépôt est bien dans `$HOME/projectSE` (chemin codé en dur) et que `./run.sh` a été lancé une fois à la main.
- **Exécution sous Windows :**  
  Le pipeline a été testé sous Ubuntu 24.04. Sous Windows, utiliser **WSL 2 (Ubuntu)**. Git Bash n'est pas garanti : le script cherche le Python du venv dans `.venv/bin/`, alors que Git Bash crée `.venv/Scripts/`.

---

## 10. Corrections issues du test croisé

Le Groupe 18 a testé le pipeline sur Ubuntu 24.04 (VMware) : démarrage, exécution complète et livrables **OK**, avis global **validé**. Leurs remarques ont été intégrées dans ce README :

| Remarque du test croisé | Correction |
| :--- | :--- |
| Les tests échouent sans l'environnement virtuel (`ModuleNotFoundError: matplotlib`) | `source .venv/bin/activate` ajouté avant les tests (section 6) et avant l'exécution modulaire (section 5) |
| `python3-venv` absent d'une Ubuntu vierge | Commande `apt install` ajoutée (section 4) et entrée FAQ |
| Chemin `projectSE` codé en dur dans le Cron | Avertissement et procédure d'adaptation ajoutés (sections 4 et 8, FAQ) |

---

## Membres de l'Équipe et Contributions

**Groupe 11 — Répartition des responsabilités :**

- **AGGOUN SOUNDOUS YASMINE — Shell :** aide et options de `run.sh`, gestion des erreurs de fichiers, commandes Shell avec pipes et filtrage des variants.
- **DAOUDI ANFAL — Programmation de script :** gestion des erreurs et des signaux avec `trap`, mode sécurisé `set -Eeuo pipefail`, fonctions réutilisables et relance du pipeline.
- **TOUBAL ABDELGHANI — Python :** analyse statistique avec `argparse`, génération des figures, export CSV/JSON et tests automatiques.
- **ZAOURAR MERIEM — Administration système :** environnement virtuel, dépendances, tâche Cron, journaux horodatés et permissions.
- **BELFERKOUS IMENE — Fichiers et données :** organisation de l'arborescence, sommes SHA-256, archivage daté des résultats et documentation des sources.
