### Installation:

1. The software requirements are managed by the **`conda`** package manager. Please install **`conda`** as described at [https://conda.io/docs/install/quick.html](https://conda.io/docs/install/quick.html). 

```
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh
```
2. Download the **nanopore-transcriptome-analysis** pipeline into a folder named `nanopore-transcriptome-analysis`.  <!--  This tutorial requires the **`git-lfs`** large file support capabilities which should be installed through **`conda`** first.     conda install -c conda-forge git-lfs --> 
``` 
    git clone https://github.com/martynakgajos/nanopore-transcriptome-analysis.git nanopore-transcriptome-analysis
```
3. Change your working directory into the new `nanopore-transcriptome-analysis` folder
```
    cd nanopore-transcriptome-analysis
```
4. Install conda software dependencies with
```
    conda env create --name nanopore-transcriptome-analysis --file environment.yaml
```
5. Initialise conda environment with 
```
    source activate nanopore-transcriptome-analysis
```

### Usage: 

In your Conda environment, and in the nanopore-transcriptome-analysis directory,

1. Edit the provided **`config.yaml`** file to match your own study design
2. Run the Snakefile workflow (the command assumes 1 available thread; adjust to match your computer's capabilities)
```
    snakemake -j 1 all
```

The pipeline is based on the **ont_tutorial_pinfish**.

