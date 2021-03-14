#!/usr/bin/env python3

"""
File name: fetch_features.py
Author: Ivan Munoz-Gutierrez
Date created: 07/26/2020
Date last modified: 03/13/2021
Python version: 3.9
"""

import csv
import sys
from Bio import Entrez
import cs50
from database import make_list_accessions, mk_uid_batch_list, parser
from database import clean_features, get_biosample_numbers


def user_input():
    """
    Gets and processes user input

    Return
    ------
    input_data : dictionary
        Information provided by user.

        Example of return dictionary:
        input_data = {
            "type_list": "accession",
            "get_biosample": True,
            "infile": "accession_list.txt",
            "email_address": "user@email.com"
        }
        In this example, the user is providing a list of accession numbers and
        wants to access the BioSample number associated to every accession
        number of the provided list.
    """
    # Dictionary to save user input
    input_data = {}
    # Checking the correct useage of the program
    if len(sys.argv) != 1:
        sys.exit("usage: python fetch_features.py")

    # Getting type of input data.
    while True:
        type_list = cs50.get_string(
            "\nDoes your list have accession numbers or "
            "biosample numbers (type accession or biosample)? ")
        type_list = type_list.lower()
        input_data["type_list"] = type_list
        if type_list == 'biosample':
            input_data["get_biosample"] = True
            break
        if type_list == 'accession':
            break

    # If list have accession numbers ask if user wants the features of all
    # related accession numbers that belong to the same BioSample number
    if type_list == 'accession':
        # Asking about getting data from all BioSample related acc numbers
        while True:
            biosample = cs50.get_string(
                "\nIf you have a list of accession numbers, do you want to "
                "get the most updated features of \nall the related accession "
                "numbers that belong to the same BioSample (type yes or no)? ")
            biosample = biosample.lower()
            if biosample == 'yes':
                input_data["get_biosample"] = True
                break
            if biosample == 'no':
                input_data["get_biosample"] = False
                break

    # Getting name of infile.txt
    input_data["infile"] = cs50.get_string(
        '\nProvide the name of your infile.txt: ')

    # IMPORTANT: always provide your email address to GenBank
    input_data["email_address"] = cs50.get_string(
        "\nProvide your email address to the NCBI: ")

    return input_data


def fetch_from_accession(results, fields, submission_list):
    """
    Fetches features from a list of accession numbers.

    Uses Entrez.post to ask for a list of accession numbers and Entrez.fetch to
    retrieve the information of the posted list of accession numbers. After
    retrieving information, the data of every accession number is parse using
    the parser function implemented in database.py.

    Parameters
    ----------
    results : file object
        Opened output file to save the fetched features
    fields : list
        List of headers of the results table. This is the list of information
        that will be fetch from every accession number.
    submission_list : list
        A list of strings containg accession numbers separated by commas. The
        number of accession numbers in every string is determined by the
        varible batch_size of the mk_acc_num_batches function. An example of a
        submission_list with a batch size of three look like the following:
        ['CP049609.1,CP028704.1,CP043542.1',
         'CP040107.1,CP041747.1,CP042638.1',
         'CP015023.1,CP049163.1,CP051714.1']
    """
    # Create DictWriter
    writer = csv.DictWriter(results, fields)
    # Declaring end
    end = 0

    # Fetching the information from GenBank by batches
    for set_number, submission in enumerate(submission_list):
        start = end
        # submission_list is a list of accession numbers separated by
        # commas. Therefore, the number of commas indicate the number of
        # accession numbers.
        batch_size = submission.count(',') + 1
        end = end + batch_size

        # Printing download batch record
        print("Going to download record %i to %i" % (start + 1, end))

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

        # Parsing the data fetched from NCBI
        records = parser(fetch_handle, set_number + 1)

        # Saving the retrived data in the csv file
        for i in range(len(records)):
            writer.writerow(records[i])

        # Closing fetch_handle
        fetch_handle.close()


def fetch_from_biosample(results, fields, list_biosamples):
    """
    Reads a list of BioSample numbers to retrieve the features of all accession
    numbers associated to every BioSample number.

    A BioSample number can have one or more associated accession numbers. For
    example, SAMN07169263 has three accession numbers (CP049611.1, CP049610.1,
    CP049609.1). The first two accession numbers belong to plasmids and the
    last one to a chromosome. This function retrieves information from all the
    accession numbers linked to a BioSample number. In the case of SAMN07169263
    this function will fetch features of the three above mentioned accession
    number.

    Parameters
    ----------
    results : file object
        Opened output file to save the fetched features
    fields : list
        List of headers of the table that will save the results. This is the
        list of information that will be fetch from every accession number
        associated to a BioSample number.
    list_biosamples: list
        List of BioSample numbers requested by the user.

    Notes
    -----
    Some BioSample numbers are associated to accession numbers that does not
    have any relevant information as contigs of few hundreds of nucleotides.
    Also, in addition to updated accession numbers, some BioSample numbers have
    outdated accession numbers. To clean all the fetched information from NCBI
    this program uses the clean_features function provided in database.py.
    """
    # Create DictWriter
    writer = csv.DictWriter(results, fields)
    lenght_acc_list = len(list_biosamples)

    # Iterating over the list of BioSample numbers (list_biosamples)
    for query, accession in enumerate(list_biosamples):
        # Number to keep track set_number of sequences (query), it is
        # important in case the connection to NCBI is interrupted so we can
        # know where to continue downloading
        set_number = query + 1

        # Printing download record
        print(
            f"Going to download record {query + 1} "
            f"of {lenght_acc_list}")
        print(f"Accessing BioSample number: {accession}")

        # Searching for the BioSample accession number. We need usehistory
        # to get the QueryKey and the WebEnv which define our history
        # session and can be used to performe searches of data.
        search_handle = Entrez.esearch(
            db="nuccore", term=accession, usehistory="y")

        # Copying information in computer memory
        search_results = Entrez.read(search_handle)

        # Closing the handle
        search_handle.close()

        # Counting the number of results (number of sequences)
        count = int(search_results["Count"])
        print(f"Number of requested sequences from BioSample: {count}")

        # Copying cookie "WebEnv" and query "QueryKey" from our history
        # session.
        # WevEnv -> Web environment string returned from a previous
        # ESearch, EPost or ELink call; QueryKey -> Integer query key
        # returned by a previous ESearch, EPost or ELink call
        webenv = search_results["WebEnv"]
        query_key = search_results["QueryKey"]

        # Number of sequences to be requested by batch.
        # A batch of 500 is the max that we can request.
        batch_size = 500

        # biosample_records will hold all the features of accession numbers
        # associated to the BioSample number that is being analyzed
        biosample_records = []

        # Fetching information from GenBank by batches
        for start in range(0, count, batch_size):
            end = min(count, start + batch_size)
            # Printing download batch record
            print(
                f"Going to download record {start + 1} to {end} "
                f"from set_number {set_number}")
            # Getting information
            # db -> database, nuccore -> nuleotide, rettype -> retrieval type,
            # retmode -> determines the format of the return output,
            # retstart -> sequential index of the first UID in the retrieved
            # set_number to be shown in the XML output, retmax -> total number
            # of UIDs from the retrieved set_number to be shown in the XML
            # output, idtype-> specifies the type of identifier to return for
            # sequence databases, acc -> accesion number
            fetch_handle = Entrez.efetch(
                db="nuccore",
                rettype="gb",
                retmode="text",
                retstart=start,
                retmax=end,
                webenv=webenv,
                query_key=query_key,
                idtype="acc"
            )

            # Parsing the data fetched from NCBI
            records = parser(fetch_handle, set_number)

            # Extending the records obtained in biosample_records for future 
            # cleaning of the data, i.e to get the most updated data
            biosample_records.extend(records)

        # Using the clean_features function to clean data and obtain the
        # most updated information.
        # clean_features() returns a list of dictionaries
        uptaded_features = clean_features(biosample_records)

        # Saving the updated retrived data in the csv file
        for _, updated in enumerate(uptaded_features):
            writer.writerow(updated)

        print(
            f"Number of sequences saved after processing: "
            f"{len(uptaded_features)}\n")

        # Closing fetch_handle
        fetch_handle.close()


def fetch_data(type_list, get_biosample, infile, email_address):
    """
    Fetchs features from a list of accession or BioSample numbers provided in a
    txt file.

    Reads a list of UIDs of a txt file and converts the UIDs into a comma-
    separated list that will be used to retrieve information via Entrez. The
    list of the input file must have a header or the first UIDs is not going to
    be analyzed.

    Parameters
    ----------
    type_list : string
        Type of UIDs provided by the user in the input file, i.e. "accession"
        or "biosample".
    get_biosample : bool
        If the user request to get features from all the accession numbers
        associated to a BioSample number get_biosample is True. Otherwise is
        False.
    infile : string
        Name of the input infile. For example, accession_list.txt.
    email_address : string
        Email address of user as requested by NCBI.
    """
    # Providing email address to NCBI
    Entrez.email = email_address
    # Reading txt infile and creating a list of accession or BioSample numbers
    list_accessions = make_list_accessions(infile)
    # Counting the number of requested accession or BioSample numbers
    count = len(list_accessions)
    print(f"\nRequested accession or BioSample numbers: {count}\n")

    # Opening our results file to write the fetched data in csv format
    with open("results.csv", "w") as results:
        # Field names or headers for the csv table in the output file
        fields = [
            "set_batch", "description", "accession", "size",
            "molecule", "mod_date", "topology", "mol_type", "organism",
            "strain", "isolation_source", "host", "plasmid", "country",
            "lat_lon", "collection_date", "note", "serovar", "collected_by",
            "genotype", "BioProject", "BioSample", "Assem_Method",
            "Gen_Coverage", "Seq_Technol"]
        # Create DictWriter
        writer = csv.DictWriter(results, fields)
        # Writing headers
        writer.writeheader()

        # Workig with a list of accession numbers
        if type_list == 'accession' and (not get_biosample):
            # Formating list of accession numbers in batches
            submission_list = mk_uid_batch_list(list_accessions)
            # Fetching features
            print("\nFetching features\n")
            fetch_from_accession(results, fields, submission_list)

        # Working with a list of accession numbers and getting all related
        # BioSample associated accession numbers and their features
        if type_list == 'accession' and get_biosample:
            # Formating list of accession numbers in batches
            submission_list = mk_uid_batch_list(list_accessions)
            # Getting BioSample numbers
            print(
                "Retrieving BioSample numbers from list of accession numbers")
            list_biosamples = get_biosample_numbers(
                submission_list, email_address)
            # Fetching features
            print("\nFetching features\n")
            fetch_from_biosample(results, fields, list_biosamples)

        # Working with a list of BioSample numbers
        if type_list == 'biosample':
            # Fetching features
            print("\nFetching features\n")
            fetch_from_biosample(results, fields, list_accessions)

    # If everything was done OK print Done and exit the program
    print("""Done!\nYou should have a results.csv file in your folder""")



def main():
    """
    Main function
    """
    # Getting user input
    input_data = user_input()
    # Fetching information
    fetch_data(
        input_data["type_list"],
        input_data["get_biosample"],
        input_data["infile"],
        input_data["email_address"])

    sys.exit(0)

if __name__ == '__main__':
    main()
