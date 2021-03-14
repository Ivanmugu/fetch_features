"""
File name: database.py
Authos: Ivan Munoz-Gutierrez
Date created: 07/25/2020
Date last modified: 03/13/2021
Python version: 3.9

This file contains functions to be used by fetch_features.py and biosample.py
"""

from Bio import Entrez
from Bio import SeqIO
from datetime import datetime
import sqlite3


def make_list_accessions(infile):
    """
    Reads a txt file that contains a list of UIDs, as accession or BioSample
    numbers, and creates a python list with these UIDs.

    Parameters
    ----------
    infile : string
        Path to txt file that contains UIDs as accession or BioSample numbers.
        Every UID must be separated by a new line character (enter). The list
        can be created in excell by saving the file as txt or using a text
        editor. The txt file must have a header or the function will not
        include the first UID.

    Returns
    -------
    list_accessions: list
        Returns a list of accession numbers. For example, a return list from a
        txt file with nine accession number would look like the following list:
        ['CP049609.1', 'CP028704.1', 'CP043542.1', 'CP040107.1', 'CP041747.1',
         'CP042638.1', 'CP015023.1', 'CP049163.1', 'CP051714.1']
    """
    # Opening infile.txt
    with open(infile, 'r') as reader:
        # Skip the header
        next(reader)
        # Creating a list of accession numbers
        list_accessions = reader.readlines()

    # Remove the '\n' character
    for i, accession_number in enumerate(list_accessions):
        list_accessions[i] = accession_number.replace('\n', '')

    return list_accessions


def mk_uid_batch_list(list_accessions):
    """
    Makes a list of batches of commma-delimited UIDs.

    Every batch created is a set of UIDs separated by commas, and every batch
    has at least the indicaded number of UIDs provided in the batch_size
    variable. Examples of UIDs are accession numbers and BioSample numbers.

    Parameters
    ----------
    list_accessions : list
        List of UIDs to be proccessed

    Returns
    -------
    submission_list : list
        List of strings containg UIDs separated by commas. An example of a
        submission_list created with accession numbers and a batch size of
        three look like the following:
        ['CP049609.1,CP028704.1,CP043542.1',
         'CP040107.1,CP041747.1,CP042638.1',
         'CP015023.1,CP049163.1,CP051714.1']
    """
    # Number of sequences or UIDs to be proccessed
    number_seq = len(list_accessions)
    # Number of sequences or UIDs to be requested by batch. Limit the batch to
    # less than 200.
    batch_size = 5
    # Declaring the list of batches of comma-delimited UIDs that will be return
    # after processing
    submission_list = []
    # Counter to access the list_accessions
    counter_accessions = 0

    # Loop to create the list of UIDs by batches
    for start in range(0, number_seq, batch_size):
        end = min(number_seq, start + batch_size)
        # This list is going to save temporarily the batch of UIDs that are
        # going to be converted into a string of comma-separated UIDs
        set_list = []
        # Making batches
        for _ in range(start, end):
            set_list.append(list_accessions[counter_accessions])
            counter_accessions += 1
        # Converting the list into string
        set_list = ','.join(set_list)
        submission_list.append(set_list)

    return submission_list


def get_biosample_numbers(submission_list, email_address):
    """
    Retrieves BioSample numbers from a list that contains batches of accession
    numbers.

    This function uses Entrez.epost to upload batches of accession numbers to
    the Entrez History server. To avoid problems with large batches of
    accession numbers, limit the number of accession number per batch to 200.

    Parameters
    ----------
    submission_list : list
        List of strings containg batches of accession numbers separated by
        commas. The following is an example of a submission_list containing
        batches of three accession numbers per element in the list:
        ["CP049609.1,CP028704.1,CP043542.1",
         "CP040107.1,CP041747.1,CP042638.1",
         "CP015023.1,CP049163.1,CP051714.1"]
    email_address : string
        User's email address is requested by NCBI

    Returns
    -------
    biosample_numbers : list
        List of BioSample numbers retrieved from Entrez by using the provided
        list of accession numbers. The following is a return list containing
        the corresponding BioSample numbers of the example provided in the
        Parameters section:
        ["SAMN07169263", "SAMN08875353", "SAMEA104140560",
         "SAMN05360217", "SAMN12302771", "SAMN12500846",
         "SAMN04202539", "SAMN14133047", "SAMN14609782"]

    Notes
    -----
    For more information about Entrez.epost read chapter 9.4 of Biopython
    Tutorial and Cookbook (Biopython 1.76) and visit:
    https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.EPost
    """
    # Creating list to save BioSample numbers
    biosample_numbers = []

    # Provide email address to GeneBank
    Entrez.email = email_address

    # Initializing last accession number
    end = 0

    # Fetching BioSample numbers from GenBank by batches
    for _, submission in enumerate(submission_list):
        start = end
        # submission_list is a list of accession numbers separated by
        # commas. Therefore, the number of commas indicate the number of
        # accession numbers.
        batch_size = submission.count(',') + 1
        end = end + batch_size

        # Printing download batch record
        print(f"Retrieving BioSample numbers from record {start + 1} to {end}")

        # Posting the submission_list.
        # Because we are requesting information from a huge list of acc
        # numbers, we have to use the ".epost" function which uploads a
        # list of UIs (acc numbers) for use in subsequent searches.
        # From .epost we can get the QueryKey and the WebEnv which define
        # our history session and can be used to performe searches of data.
        posting = Entrez.epost('nuccore', id=submission)
        search_results = Entrez.read(posting)

        # Copying cookie "WebEnv" and query "QueryKey" from our history
        # session to keep track of our batch fetching.
        # WevEnv -> Web environment string returned from a previous
        # ESearch, EPost or ELink call; QueryKey -> Integer query key
        # returned by a previous ESearch, EPost or ELink call
        webenv = search_results["WebEnv"]
        query_key = search_results["QueryKey"]

        # Getting the batch information
        # db -> database, nuccore -> nuleotide, rettype -> retrieval type
        # retmode -> determines the format of the return output
        # retstart -> sequential index of the first UID in the retrieved
        # set_number to be shown in the XML output
        # retmax -> total number of UIDs from the retrieved set_number to be
        # shown in the XML output
        # idtype-> specifies the type of identifier to return for sequence
        # databases, acc -> accesion number
        fetch_handle = Entrez.efetch(
            db="nuccore",
            rettype="gb",
            retmode="text",
            retstart=0,
            retmax=batch_size,
            webenv=webenv,
            query_key=query_key,
            idtype="acc"
        )

        # Parsing throw the fetched information
        for seq_record in SeqIO.parse(fetch_handle, "gb"):
            # Looping over database cross-references (dbxrefs)
            for xref in enumerate(seq_record.dbxrefs):
                list_dbxrefs = xref[1].split(':')
                # Getting the BioSample number
                if 'BioSample' in list_dbxrefs:
                    biosample_numbers.append(list_dbxrefs[1])
                    break

    return biosample_numbers


def parser(fetch_handle, set_number):
    """
    Function to parse information fetched from nuccore NCBI

    Parameters
    ----------
    fetch_handle : network connection
        Pointer to a network connection that will allow to download and parse
        the requested sequences from the internet
    set_number : inteter
        Number of the batch of sequences that is being analyzed. This number
        will help to keep track of the downloaded sequences.
    seq_counter : integer
        Number of the position of the first accession number in the set. If
        set_number is 1, seq_counter has to be 1.
    
    Returns
    -------
    records : list
        List of dictionaries. Every dictionary corresponds to fetch features of
        one accession number.
    """
    # Creating a list to save records in dictionaries
    records = []
    info = {}

    # Parsing throw the fetched information
    for seq_record in SeqIO.parse(fetch_handle, "gb"):
        # keeping track of the set_number (or batch) analyzed
        info['set_batch'] = set_number

        # Extracting description
        info['description'] = seq_record.description

        # Extracting sequence id, i.e. accession number
        info['accession'] = seq_record.id

        # '.seq' is an object with the sequence itself
        info['size'] = len(seq_record.seq)

        #######################################################################
        # '.annotations' is a dictionary of aditional information about the
        # sequence as last modification date, topology, sequence_version,
        # organims, references, etc.
        #######################################################################
        # Getting date and give format
        if 'date' in seq_record.annotations:
            mod_date = seq_record.annotations['date']
            mod_date = datetime.strptime(mod_date, '%d-%b-%Y')
            mod_date = datetime.strftime(mod_date, "%Y-%m-%d")
            info['mod_date'] = mod_date
        else:
            info['mod_date'] = 'missing'
        # Getting topology
        if "topology" in seq_record.annotations:
            info['topology'] = seq_record.annotations["topology"]
        else:
            info['topology'] = "missing"

        # Checking whether it is chromosome or plasmid. The 4,000,000 used in
        # size is an arbitrary number that considers chormosomes of that length
        # Getting chromosome
        if 'chromosome' in info['description']:
            info['molecule'] = 'chromosome'
        # Getting chromosome
        elif info['size'] >= 4000000 and info['topology'] == 'circular':
            info['molecule'] = 'chromosome'
        # Getting plasmid
        elif 'plasmid' in info['description']:
            info['molecule'] = 'plasmid'
        else:
            info['molecule'] = 'missing'

        # '.features' is a list of SeqFeatures objects with more structured
        # information about the features on a sequence
        feature = seq_record.features

        # Looping throw list feature
        for index in feature:
            # '.type' is only a description of the type of feature
            # that could be source, CDS, gene, etc.
            # In source we can find organism, strain, host, country, etc.
            if index.type == "source":
                # Creating a dictionary of the qualifiers from source
                dictionary = dict(index.qualifiers)
                # '.get' gives a list
                # Getting molecule type
                if "mol_type" in dictionary:
                    info['mol_type'] = dictionary.get("mol_type")[0]
                else:
                    info['mol_type'] = 'missing'
                # Getting organism
                if 'organism' in dictionary:
                    info['organism'] = dictionary.get('organism')[0]
                else:
                    info['organism'] = 'missing'
                # Getting strain
                if 'strain' in dictionary:
                    info['strain'] = dictionary.get('strain')[0]
                else:
                    info['strain'] = 'missing'
                # Getting isolation source
                if 'isolation_source' in dictionary:
                    info['isolation_source'] =\
                        dictionary.get('isolation_source')[0]
                else:
                    info['isolation_source'] = 'missing'
                # Getting host
                if 'host' in dictionary:
                    info['host'] = dictionary.get('host')[0]
                else:
                    info['host'] = 'missing'
                # Getting plasmid
                if 'plasmid' in dictionary:
                    info['plasmid'] = dictionary.get('plasmid')[0]
                else:
                    info['plasmid'] = 'missing'
                # Getting country
                if 'country' in dictionary:
                    info['country'] = dictionary.get('country')[0]
                else:
                    info['country'] = 'missing'
                # Getting coordinates
                if "lat_lon" in dictionary:
                    info['lat_lon'] = dictionary.get("lat_lon")[0]
                else:
                    info['lat_lon'] = 'missing'
                # Getting collection date
                if "collection_date" in dictionary:
                    info['collection_date'] =\
                        dictionary.get("collection_date")[0]
                else:
                    info['collection_date'] = 'missing'
                # Getting notes
                if "note" in dictionary:
                    info['note'] = dictionary.get("note")[0]
                else:
                    info['note'] = 'missing'
                # Getting serovar
                if "serovar" in dictionary:
                    info['serovar'] = dictionary.get("serovar")[0]
                else:
                    info['serovar'] = 'missing'
                # Getting collected by
                if "collected_by" in dictionary:
                    info['collected_by'] = dictionary.get("collected_by")[0]
                else:
                    info['collected_by'] = 'missing'
                # Getting genotype
                if "genotype" in dictionary:
                    info['genotype'] = dictionary.get("genotype")[0]
                else:
                    info['genotype'] = 'missing'

                break

        # '.dbxrefs' is a list populated from any PROJECT or DBLINK
        # Checking if .dbxrefs has any information
        if len(seq_record.dbxrefs) == 0:
            info['BioProject'] = 'missing'
            info['BioSample'] = 'missing'

        # Converting the list dbxrefs into dictionary
        dictionary_dbxrefs = {}
        for i in range(len(seq_record.dbxrefs)):
            s = seq_record.dbxrefs[i].split(":")
            dictionary_dbxrefs[s[0]] = s[1]

        # Saving BioProject and BioSample
        if "BioProject" in dictionary_dbxrefs:
            info['BioProject'] = dictionary_dbxrefs.get("BioProject")
        else:
            info['BioProject'] = 'missing'
        if "BioSample" in dictionary_dbxrefs:
            info['BioSample'] = dictionary_dbxrefs.get("BioSample")
        else:
            info['BioSample'] = 'missing'

        # Getting the Genome-Assembly-Data
        # Checking if the sequence has structured_comment
        if (("structured_comment" in seq_record.annotations)
            and ("Genome-Assembly-Data" in (
                seq_record.annotations["structured_comment"]))):
            # Getting Assembly Method
            if ("Assembly Method" in (
                    seq_record.annotations["structured_comment"]
                                          ["Genome-Assembly-Data"])):
                info['Assem_Method'] = (
                    seq_record.annotations["structured_comment"]
                                          ["Genome-Assembly-Data"]
                                          ["Assembly Method"])
            else:
                info['Assem_Method'] = "missing"
            # Getting Genome Coverage
            if "Genome Coverage" in (
                seq_record.annotations["structured_comment"]
                                      ["Genome-Assembly-Data"]):
                info['Gen_Coverage'] = (
                    seq_record.annotations["structured_comment"]
                                          ["Genome-Assembly-Data"]
                                          ["Genome Coverage"])
            else:
                info['Gen_Coverage'] = "missing"
            # Getting Sequencing Technology
            if "Sequencing Technology" in (
                seq_record.annotations["structured_comment"]
                                      ["Genome-Assembly-Data"]):
                info['Seq_Technol'] = (
                    seq_record.annotations["structured_comment"]
                                          ["Genome-Assembly-Data"]
                                          ["Sequencing Technology"])
            else:
                info['Seq_Technol'] = "missing"
        # Extracting data from structured comment
        elif "structured_comment" in seq_record.annotations and (
            "Assembly-Data" in (
                seq_record.annotations["structured_comment"])):
            # Getting Assembly Method
            if "Assembly Method" in (
                seq_record.annotations["structured_comment"]
                                      ["Assembly-Data"]):
                info['Assem_Method'] = (
                    seq_record.annotations["structured_comment"]
                                          ["Assembly-Data"]
                                          ["Assembly Method"])
            else:
                info['Assem_Method'] = "missing"
            # Getting Genome Coverage
            if "Genome Coverage" in (
                seq_record.annotations["structured_comment"]
                                      ["Assembly-Data"]):
                info['Gen_Coverage'] = (
                    seq_record.annotations["structured_comment"]
                                          ["Assembly-Data"]
                                          ["Genome Coverage"])
            else:
                info['Gen_Coverage'] = "missing"
            # Getting Sequencing Technology
            if "Sequencing Technology" in (
                seq_record.annotations["structured_comment"]
                                      ["Assembly-Data"]):
                info['Seq_Technol'] = (
                    seq_record.annotations["structured_comment"]
                                          ["Assembly-Data"]
                                          ["Sequencing Technology"])
            else:
                info['Seq_Technol'] = "missing"

        else:
            info['Assem_Method'] = "missing"
            info['Gen_Coverage'] = "missing"
            info['Seq_Technol'] = "missing"

        # Copying info (dictionary) into records (list)
        records.append(info.copy())

    return records


def clean_features(raw_data):
    """
    Cleans up the features obtained from a BioSample number and get the most
    updated information

    Parameters
    ----------
    raw_data : list
        List of dictionaries holding all the data downloaded from a specific
        BioSample associated accession numbers.

    Returns
    -------
    recods : list
        List of dictionaries holding the most updated information from
        BioSample associated accession numbers.
    """

    # Creating empty features.db
    open('features.db', 'w').close()

    # Opening features.db for sqlite3
    conn = sqlite3.connect('features.db')

    # Creating a cursor
    c = conn.cursor()

    # Creating table in SQL for the features of the raw data fetched from NCBI
    c.execute("""CREATE TABLE features_raw (
                set_batch int,
                description text,
                accession text,
                size integer,
                molecule text,
                mod_date text,
                topology text,
                mol_type text,
                organism text,
                strain text,
                isolation_source text,
                host text,
                plasmid text,
                country text,
                lat_lon text,
                collection_date text,
                note text,
                serovar text,
                collected_by text,
                genotype text,
                BioProject text,
                BioSample text,
                Assem_Method text,
                Gen_Coverage text,
                Seq_Technol text
                )""")
    conn.commit()

    # Transfering results obtained from NCBI into the table features_raw
    # created in SQL
    for _, raw_result in enumerate(raw_data):
        c.execute("""INSERT INTO features_raw
                     VALUES(:set_batch, :description, :accession,
                            :size, :molecule, :mod_date, :topology, :mol_type,
                            :organism, :strain, :isolation_source, :host,
                            :plasmid, :country, :lat_lon, :collection_date,
                            :note, :serovar, :collected_by, :genotype,
                            :BioProject, :BioSample, :Assem_Method,
                            :Gen_Coverage, :Seq_Technol)""",
                  {'set_batch': int(raw_result['set_batch']),
                   'description': raw_result['description'],
                   'accession': raw_result['accession'],
                   'size': int(raw_result['size']),
                   'molecule': raw_result['molecule'],
                   'mod_date': raw_result['mod_date'],
                   'topology': raw_result['topology'],
                   'mol_type': raw_result['mol_type'],
                   'organism': raw_result['organism'],
                   'strain': raw_result['strain'],
                   'isolation_source': raw_result['isolation_source'],
                   'host': raw_result['host'],
                   'plasmid': raw_result['plasmid'],
                   'country': raw_result['country'],
                   'lat_lon': raw_result['lat_lon'],
                   'collection_date': raw_result['collection_date'],
                   'note': raw_result['note'],
                   'serovar': raw_result['serovar'],
                   'collected_by': raw_result['collected_by'],
                   'genotype': raw_result['genotype'],
                   'BioProject': raw_result['BioProject'],
                   'BioSample': raw_result['BioSample'],
                   'Assem_Method': raw_result['Assem_Method'],
                   'Gen_Coverage': raw_result['Gen_Coverage'],
                   'Seq_Technol': raw_result['Seq_Technol']})
    conn.commit()

    # Getting the most updated files and ordering by size
    c.execute(
        """
        SELECT *
          FROM (
               SELECT table_one.*
                 FROM features_raw table_one
                WHERE table_one.mod_date = (
                      SELECT MAX(table_two.mod_date)
                        FROM features_raw table_two
                       WHERE table_two.size = table_one.size))
         WHERE molecule != "missing"
         ORDER BY size DESC
        """)

    # Extracting the updated results from SQL. The fetchall will create a list
    # of tuples
    results = c.fetchall()

    # Converting results (list of tuples) into a list of dictionaries
    records = []
    info = {}
    for _, updated_result in enumerate(results):
        info['set_batch'] = int(updated_result[0])
        info['description'] = updated_result[1]
        info['accession'] = updated_result[2]
        info['size'] = int(updated_result[3])
        info['molecule'] = updated_result[4]
        info['mod_date'] = updated_result[5]
        info['topology'] = updated_result[6]
        info['mol_type'] = updated_result[7]
        info['organism'] = updated_result[8]
        info['strain'] = updated_result[9]
        info['isolation_source'] = updated_result[10]
        info['host'] = updated_result[11]
        info['plasmid'] = updated_result[12]
        info['country'] = updated_result[13]
        info['lat_lon'] = updated_result[14]
        info['collection_date'] = updated_result[15]
        info['note'] = updated_result[16]
        info['serovar'] = updated_result[17]
        info['collected_by'] = updated_result[18]
        info['genotype'] = updated_result[19]
        info['BioProject'] = updated_result[20]
        info['BioSample'] = updated_result[21]
        info['Assem_Method'] = updated_result[22]
        info['Gen_Coverage'] = updated_result[23]
        info['Seq_Technol'] = updated_result[24]
        records.append(info.copy())
    conn.close()

    return records
