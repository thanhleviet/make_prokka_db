# Download and Make a database for use with Prokka

The script `make_prokka_db.py` can be used to:

    - Download complete genomes in genbank format from NCBI given a list of NCBI taxonomy IDs.
    - Build a blast database for use with Prokka.

```bash
./make_prokka_db.py --help
usage: make_prokka_db.py [-h] [-t TAXID] [-n NAME] [-o OUTDIR]
                         [-s {refseq,genbank}]
                         [-l {all,complete,chromosome,scaffold,contig}]
                         [-g GROUP] [-p N] [-e EXT] [-b]

Download and Make a database for use with Prokka

optional arguments:
  -h, --help            show this help message and exit
  -t TAXID, --taxid TAXID
                        Only download sequences of the provided NCBI taxonomy
                        ID. A comma-separated list of taxids is also possible.
                        For example: "9606,9685". (default: 93071)
  -n NAME, --name NAME  A name for the database (default: new_prokka_database)
  -o OUTDIR, --outdir OUTDIR
                        A directory for storing intermediate outputs (default:
                        /home/ubuntu/mydata/salmonella/prokka_db_maker/ncbi)
  -s {refseq,genbank}, --section {refseq,genbank}
                        NCBI section to download (default: genbank)
  -l {all,complete,chromosome,scaffold,contig}, --assembly-level {all,complete,chromosome,scaffold,contig}
                        Assembly level of genomes to download (default:
                        complete)
  -g GROUP, --group GROUP
                        Taxonomic group, i.e bacteria, viral, etc (default:
                        bacteria)
  -p N, --parallel N    Run N downloads and converting gbk to faa in parallel
                        (default: 1)
  -e EXT, --ext EXT     File extension for scanning with sequence folder
                        (default:gz) (default: gz)
  -b, --build           Build database given from a folder of complete genbank
                        files? (default: False)
```

## Dependencies

Before running the script, please make sure you have the following dependencies: cd-hit, blast, ncbi-genome-download, biopython

External dependencies can be installed from bioconda

```bash
conda install -c conda-forge -c bioconda cd-hit blast
```

Python dependecies can be install via pip

```bash
pip install ncbi-genome-download==0.2.8 biopython
```

## Example

### Download complete sequences for *Salmonella enterica subsp. enterica serovar Typhimurium*

Taxonomy ID: 90371

```bash
./make_prokka_db.py -t 90371 -n salmonella_90371 -p 4
```

Results if run sucessfully.

```bash
Start downloading 90371
Location /home/ubuntu/mydata/salmonella/prokka_db_maker/ncbi
Start building DB for: 90371
Database salmonella_90371 for use with prokka has been saved to /home/ubuntu/mydata/salmonella/prokka_db_maker/salmonella_90371
Finished!
```

Output files:

```bash
salmonella_90371
├── salmonella_90371.faa
├── salmonella_90371.phr
├── salmonella_90371.pin
└── salmonella_90371.psq
salmonella_90371.meta
```

### Download complete sequences for *Salmonella enterica subsp. enterica serovar Typhi* and *Salmonella enterica subsp. enterica serovar Typhimurium*

Taxonomy ID: 590

```bash
./make_prokka_db.py -t 90370,90371 -n salmonella -p 4
```

### Build database from a folder (e.g. ncbi) of genbank files

```bash
./make_prokka_db.py -o ncbi -n salmonella -p 4  -b
```
