LOFAR EoR Demonstration pipeline
========================================================

lofar-eor-demo is a demonstration pipeline for the Epoch of Reionization experiment with [Lofar-EoR](http://www.lofar.org/astronomy/eor-ksp/epoch-reionization.html). It utilizes tools developed by multiple members of the LOFAR collaboration to perform data reduction from pre-processing to production of final power spectra. The main tools include:

- [DP3](https://github.com/lofar-astron/DP3)
- [AOflagger](https://sourceforge.net/p/aoflagger/wiki/Home/)
- [WSClean](https://sourceforge.net/p/wsclean/wiki/Home/)
- [Pspipe](https://gitlab.com/flomertens/pspipe)

Dependencies
------------

- [Nextflow](https://www.nextflow.io/docs/latest/index.html)
- A singularity container with all dependencies


Downloads
--------------
Download a LOFAR-EoR test observation and the singularity image from our `astro.rug.nl` cloud. This will download a `codex.zip` file into your chosen directory `/home/myname/mydir/testdata`. The `codex.zip`. Unzip the file, which will give you 2 files: the container image `leor_tools_050625.sif`, and the lofar test data `L192832.tar`

```
# Get the container file and the test data
wget --content-disposition "https://cloud.astro.rug.nl/index.php/s/ozZ3ZEaPmdTkR9y/download"
# Unzip the folder in your chosen directoy:
unzip codex.zip
```

Documentation
--------------

Once you have nextflow installed and made it executable, you can run the LOFAR-EoR tutorial by providing the path to your `.sif` image as configurable paramter in the `nextflow run` command.

```
 nextflow run /path/to/repo/lofar-eor-demo/main.nf --data.ms_files.raw '/home/myname/mydir/testdata/L192832_SAP000_SB???_uv_3c196.MS.5ch8s' --data.path '/home/myname/mydir/testdata/codex/' --average.lta_to_di.timestep 1 --average.lta_to_di.freqstep 1 --container_filename /home/myname/mydir/testdata/codex/leor_tools_050625.sif
```
