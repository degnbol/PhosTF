This is the distribution of NetworKIN 3.0. November, 2013

Basic information:

The distribution contains the ANSI-C source code of NetPhorest, a python script of NetworKIN 3.0, and several data files. For more information of this program, see kinomexplorer.info hosted by the LindingLab at CBS/DTU.


Installation:

	1) Compile netphorest
		cc -O3 -o netphorest netphorest.c -lm
	2) Install blast (ftp://ftp.ncbi.nih.gov/blast/executables/release/2.2.17/)
	NetworKIN currently supports blast version 2.2.17 or earlier. To use later versions which was named blast plus, NetworKIN.py needs to be modified accoring to the changes in blast plus.
	3) Install python (version 2.6 or 2.7)
	4) Edit the first line of NetworKIN.py as the path to the python executable installed
		EX) #!/usr/bin/env python


Running NetworKIN:

	1) open Terminal
	2) cd NetworKIN directory
	3) Run NetworKIN.py 
		Usage: ./NetworKIN.py -n path_to_netphorest_executable -b path_to_blastall Organism FastaFile SitesFile
		e.g.) ./NetworKIN.py -n /Works/NetPhorest/netphorest -b /Packages/blast-2.2.17/bin/blastall 9606 substrate_proteins.fas phospho_residues.res
	NOTE: You can store the output to a file by using UNIX output redirection (‘>’).
		.networkin.py Organism FastaFile SitesFile > OutputFile

The NetworKIN algorithm takes taxon code as organism (human: 9606, yeast: 4932). NetworKIN currently supports human and yeast, but can also be deployed on mouse data. To use NetworKIN for other organisms, see the help page of the web interface (kinomexplorer.info).


Input format:

	* FastaFile: This file is a standard fasta format file (for more information about fasta format, see http://en.wikipedia.org/wiki/FASTA_format), which contains protein names and sequences. Note that proteins names in this fasta file should correspond to protein names in the SitesFile.
	* SitesFile: NetworKIN takes three formats of SiteFile, which is automatically detected in the script.
	i) NetworKIN site file: This is a simple tab-delimited text file that contains phosphorylated residue information as the following format.
	protein name(tab)position(tab)amino acid
	ii) Direct output file from MaxQuant
	iii) Protein IDs and peptides information extracted from the output file of ProteomeDiscoverer.
	Protein names should correspond to the names in the FastaFile, and positions and amino acid are positions and amino acids of a phosphorylated residue in the corresponding sequence in the FastaFile.


Output format:

In the output of NetworKIN, each line contains a prediction result between a pair of a phosphorylation site and a kinase (KIN)/phosphatase (PTP)/phospho-binding domain.

1) Name: The name of the target protein provided by the user in the FastaFile and SiteFile as well
2) Position: The position of phosphorylation site in the sequence.
3) Tree: A tree name that the enzyme protein belongs to
4) NetPhorest Group: A name of protein group that the enzyme belongs to. The NetPhorest groups were determined by NetPhorest training pipeline.
5) Kinase/Phosphatase/Phospho-binding domain: A name of the enzyme of the prediction. This protein can be a kinase, phosphatase, or other phospho-binding domain proteins.
6) NetworKIN score: An integrated score of NetPhorest probability and STRING score
7) NetPhorest probability: Posterior probability from NetPhorest 
8) STRING score: A score that represents proximity of the target (substrate) and the enzyme in the STRING network, which is calculated by the best path search algorithm and the hub/length penalty scheme 
9) Target STRING ID: STRING ID of the target
10) Kinase/Phosphatase/Phospho-binding domain STRING ID: STRING ID of the enzyme
11) Target description: A description of the target protein
12) Kinase/Phosphatase/Phospho-binding domain description: A description of the enzyme
13) Target Name: A conventional name of the target protein
14) Kinase/Phosphatase/Phospho-binding domain Name: A conventional name of the enzyme
15) Peptide sequence window: The peptide sequence surrounding the phosphorylation site that was used to calculate NetPhorest probability
16) Intermediate nodes: A list of proteins that the best path passes through in STRING network


