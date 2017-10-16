# annopipe - Automated pipeline for annotation and parsing of predicted genes files.

This pipeline is interesting for annotation and parsing of large datasets, using this in protein sequence files obtained after the gene prediction.

--------------------------------------------------------------------------

**This version accounts resources to allow annotation with:**

[KEGG](http://kegg.jp/)

[COG](https://www.ncbi.nlm.nih.gov/COG/)

[Uniprot](http://www.uniprot.org/)

[EGGNOG](http://eggnogdb.embl.de/#/app/home)

[PFAM](http://pfam.xfam.org/)

[CAMERA](http://camera.calit2.net/)

---------------------------------------------------------------------------

1. Before Installing:

Make sure that you have installed: [HMMER](http://hmmer.org/); [Diamond](https://github.com/bbuchfink/diamond); Python && Biopython.

Before starting you also should download the last update of the databases above mentioned to a folder of your preference.

Then generates the diamond databases using the following command:

```
$ diamond makedb --in <database.fasta> -d <database.output>
```

2. Installation:

After download the files, you should update the paths in the script file. Take attention to:

A. Substitute the paths of the databases pre-formatted in diamond database;

B. HMMER databases are ready to use, just update the paths

C. **Download the program libraries from zenodo repository: [Link](https://zenodo.org/record/1013521#.WeTepCdLfCI)**
 
D. type this command:

```
$ tar -zxvf /path/to/ref_folder.tar.gz
```

E. add the ref_folder path to the script

F. add the pathway of anno_pipe.sh to your main PATH


3. Running:

```
$ anno_pipe.sh <input_fasta_file> <output_folder> <threads>
```

    H++ :: In case of threads number not be specified then the anno_pipe will use 4 threads.

--------------------------------------------------------------------------

**Specifications**

System : Linux

Uses cutoff of 50% coverage, minimum score of 50, 45% of identity and 1e-5 in Diamond-Blastp approaches, while uses 60% coverage and 1e-5, with no domains overlapping for HMM profiles.

For more details look directly the sections parsing in the script. The original searching files are kept and allows a new parsing as desired by the user.

--------------------------------------------------------------------------

**Author**: Célio Dias Santos Júnior

*Biotechonologist; Master in Molecular Biology and Genetics; PhD student in Bioinformatics*

[**Contact**](celio.diasjunior@gmail.com)

--------------------------------------------------------------------------

## License

This software is released under GNU license, also following this distribution. 
