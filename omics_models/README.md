# Omics Models

### Overview

This README describes the code used in Figure 5 of the article.

### Prerequisites

The following dependencies were used to run AlphaMissense. Any Python packages can be installed using the [requirements.txt](./requirements.txt)

* aria2
* HMMER
* Python 3.x

Download Uniprot FASTA (canonical and non-canonical) sequences for homo sapiens into tsv format if you need the dictionary produced by [00_create_uniprot_dict.py](./00_create_uniprot_dict.py)

Download desired databases to run AMS as described by google [here](https://github.com/google-deepmind/alphafold/tree/main/afdb)

### Installation

Clone the AlphaMissense [repository](https://github.com/google-deepmind/alphamissense/tree/main) to your desired path.

```
git clone https://github.com/google-deepmind/alphamissense.git
cd ./alphamissense
```

Create a Python virtual environment to run AlphaMissense and other scripts. From the AMS repository cloned in the prior step, run the following:

```
python3 -m venv ./venv
venv/bin/pip install -r requirements.txt
venv/bin/pip install -e .
```

Move [01_run_AM.py](./01_run_AM.py) to the ```alphamissense/test``` directory to run the variant of AMS.

In ```alphamissense/alphamissense/model/config.py```, set ```zero_init``` to ```False```.

Changing this parameter initializes the weights to linear instead of being zeroed out for Haiku.

Move the AMS repository cloned/downloaded to your desired directory and training and test variant text files.

```
    mv training_variants.txt <path to omics project>/omics_models/train/training_variants.txt
    mv testing_variants.txt <path to omics project>/omics_models/test/testing_variants.txt
```

Fill out desired parameters and use command as shown in [01_run_AM_batch.sh](./01_run_AM_batch.sh)

### Auxiliary Scripts

The auxiliary scripts used to create the other panels in the figure are provided for reference and should have absolute/relative paths replaced according to your local machine. These scripts utilize the numpy arrays produced by the AlphaMissense script provided in this repository.