# fetch_features

In order to run this script you need python3.
Additionaly, you need the Biopython and cs50 modules.
Therefore, you have to install the Biopython and cs50 modules in your computer.
If you want to know if you have those modules type the following in the terminal:

```bash
pip3 list
```

If you don't have those modules, you can install them by typing:
```bash
pip3 install biopython
pip3 install cs50
```

You can get more information about those modules by visiting:
* https://biopython.org/
* https://github.com/cs50/python-cs50

Finally, you need the script database.py (that is provided in this project) in the
same folder where you will run fetch_features.py

# About fetch_features.py
## Function:
Fetch information from a list of Genebank accession numbers or a list of BioSample
numbers.

The program reads the accession or BioSample numbers of the list. Then, it enters
into the nuccore database and collects all the molecular features of each accession
or BioSample number in the list.

You can create one list of accession or BioSample numbers in Excel by saving the
file as txt. The list needs a header, if it doesn't have one the first accession or
BioSample numbers is not going to be included.

As mentioned before, the program can use accession of BioSample numbers. Therefore,
when you run the program, the program will ask you:
```
Does your list have accession numbers or biosample numbers (type accession or biosample)? 
```

If you choose biosample, you will get the features of all the accession numbers
associated with the BioSample number. For example, if your list has the next two
BioSample numbers, SAMN07169263 and SAMN08875353, you will get the features of the 
following accession numbers:
```
SAMN07169263 will return the features of: CP049609.1, CP049611.1 and CP049610.1.
SAMN08875353 will return the features of: CP028704.1 and CP028705.1.
```

If you choose accession, the program will ask you to select between two options as follows:

```
If you have a list of accession numbers, do you want to get the most updated features of 
all the related accession numbers that belong to the same BioSample (yes or no)?
```

If you answer no, you will get only the features of the provided accession list. For example,
if your list has the next two accession numbers, CP049609.1 and CP028704.1, you will get the
features only of those two accession numbers.

If you answer yes, the program will first get the BioSample number of every provided accession
number in the list. Then, these BioSample numbers will give access to all the accesion numbers
associated with every BioSample number. Finally, the program will fetch the molecular features
of all the accession numbers associated with every BioSample number. For example, if your list
has the next two accession numbers, CP049609.1 and CP028704.1, the program will follow the
next steps:
```
1. Get the corresponding BioSample numbers of CP049609.1 and CP028704.1 which are SAMN07169263
and SAMN08875353, respectively.
2. Access each BioSample number, get all the associated accession numbers to each BioSample
number and get the features of all the accession numbers. Therefore, the program will get the
features of the following accession numbers:
    i)  CP049609.1, CP049611.1 and CP049610.1 that are associated to SAMN07169263.
    ii) CP028704.1 and CP028705.1 that are associated to SAMN08875353.
```

Sometimes, one BioSample number can have redundant information. For example, a plasmid can have
two accession numbers, one accesssion could have the most updated information while the second
one could have outdated information. Additionaly, the same BioSample number can have accession
numbers that carry no relevant information as short contigs of a few nucleotides. To get only
relevant information, the program removes outdated and non-relevant information.

## Example of usage
Before running the program you need a folder with all the relevant information. For example:
```
Documents
+-------- Results
          +------ fetch_features.py
          +------ database.py
          +------ accession_numbers_list.txt
```

The accession_numbers_list.txt has to have a header. The file should look like this:
```
Accession_number
CP043542.1
CP049163.1
```

Next, you run the code as follows:
```bash
python3 fetch_features.py
```

As mentioned before, after running the code you will be asked to provide some information.
In this example, we will request to get the features of all accession numbers associated to
their corresponding BioSample number. The last question is your email address. NCBI needs
your email address.
```
Does your list have accession numbers or biosample numbers (type accession or biosample)? accession

If you have a list of accession numbers, do you want to get the most updated features of 
all the related accession numbers that belong to the same BioSample (yes or no)? yes

Provide the name of your infile.txt: accession_list.txt

Provide your email address to the NCBI: user@email.com
```

After connecting to NCBI you will see the followin messages during the fetching proccess.
```
Requested accession or BioSample numbers: 2

Retrieving BioSample numbers from list of accession numbers
Retrieving BioSample numbers from record 1 to 2

Fetching features

Going to download record 1 of 2
Accessing BioSample number: SAMEA104140560
Number of requested sequences from BioSample: 204
Going to download record 1 to 204 from set_number 1
Number of sequences saved after processing: 4

Going to download record 2 of 2
Accessing BioSample number: SAMN14133047
Number of requested sequences from BioSample: 616
Going to download record 1 to 500 from set_number 2
Going to download record 501 to 616 from set_number 2
Number of sequences saved after processing: 2

Done!
You should have a results.csv file in your folder
```

That's it!
If you opened the results.csv file in Excell you will see something like this:

![alt text](https://github.com/Ivanmugu/fetch_features/blob/main/results.png)
