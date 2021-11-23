This is the distribution of NetPhorest yeast v2.1, July 2017

The distribution contains the ANSI-C source code of NetPhorest. For more information about this program, see http://kinomexplorer.info hosted by the LindingLab.

If you use this program in your research, please cite the following article: 
	* Horn et al., KinomeXplorer: an integrated platform for kinome biology studies. Nature Methods 2014 Jun;11(6):603–4. <https://doi.org/10.1038/nmeth.2968>


Installation:
	1) unzip binary_yeast_release.zip -d binary_yeast_release
	2) cd binary_human_release
	3) Compile netphorest:
		cc -O3 -o netphorest netphorest.c -lm


Usage: cat FastaFile | ./netphorest

	FastaFile: sequences in the standard fasta format (for more information about fasta format, see http://en.wikipedia.org/wiki/FASTA_format), which contains protein names and sequences. Note that proteins names in this fasta file should correspond to protein names in the SitesFile.

        NOTE: You can store the output to a file by using UNIX output redirection (‘>’).
                Usage: cat FastaFile | ./netphorest > OutputFile

