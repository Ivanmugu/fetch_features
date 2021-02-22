# fetch_features

In order to run this script you need python3.
Additionaly, you need the Biopython and cs50 modules.
Therefore, you have to install the Biopython and cs50 modules in your computer.
If you want to know if you have those modules type the following in the terminal:
pip3 list

If you don't have those modules, you can install them by typing:
pip3 install biopython
pip3 install cs50

You can get more information about those modules by visiting:
https://biopython.org/
https://github.com/cs50/python-cs50

# About fetch_features.py
Function:
Fetch information from a list of Genebank accession numbers or a list of BioSample
numbers.

The program reads the accession or BioSample numbers of the list. Then, it enters
into the nuccore database and collects all the molecular features of each accession
of BioSample number in the list.

If you provide a list of accession numbers you have two options. The first
option is to get only the features of the provided accession list. In the second
option, the program gets the BioSample number to the provided accession number.
Then, this BioSample number gives access to all the accesion numbers associated
with the BioSample number. Finally, the program fetches the molecular features of
all the accessin numbers associated with this BioSample number. In this second
option, the program selects the most updated information of every molecule
(chromosome and/or plasmid(s)).

You can create one list of accession or BioSample numbers in Excel by saving the
file as txt. The list needs a header, if it doesn't have one the first accession or
BioSample numbers is not going to be included.

# Example of usage

Usage: python3 fetch_features.py
Does your list have accession numbers or biosample numbers (accession or biosample)? accession
If you have a list of accession numbers, do you want to get the most updated features of 
all the related accession numbers that belong to the same BioSample (yes or no)? yes
Provide the name of your infile.txt: accession_list.txt
Provide your email address to the NCBI: user@email.com
Number of requested sequences: 2

Requested accesion number: CP043542.1
Corresponding BioSample number: SAMEA104140560
Number of requested sequences from BioSample: 204
Going to download record 1 to 204 from set 1
Number of sequences saved after processing: 4

Requested accesion number: CP049163.1
Corresponding BioSample number: SAMN14133047
Number of requested sequences from BioSample: 616
Going to download record 1 to 500 from set 2
Number of sequences saved after processing: 0
Going to download record 501 to 616 from set 2
Number of sequences saved after processing: 2
Done!
You should have a results.csv file in your folder
