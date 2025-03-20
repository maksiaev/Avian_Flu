
import os
import glob 
import pandas as pd

# Function to convert fasta file to dataframe 
def fasta_df(file_name):

    fasta = pd.DataFrame()
    headers = []
    isolate_ids = []
    isolate_names = []
    subtypes = []
    segments = []
    collection_dates = []
    sequences = []
    with open(file_name) as f:
        lines = f.readlines()
        for num, line in enumerate(lines):
            if line[0] != ">": # If it's not a header
                if lines[num - 1].strip() not in headers: # And it's not a header we've seen before
                    header = lines[num - 1][1:].strip() # Remove the ">"
                    headers.append(header) 
                    isolate_ids.append(header.split("|")[0])
                    isolate_names.append(header.split("|")[1]) # We'll need to extract data from this too
                    subtypes.append(header.split("|")[2])
                    segments.append(header.split("|")[3])
                    if header.split("|")[4] == "2024-01-01":
                        collection_dates.append("2024") # No samples were collected 1/1/2024, these are all unknown 
                    elif header.split("|")[4] == "2025-01-01":
                        collection_dates.append("2025")
                    else: 
                        collection_dates.append(header.split("|")[4])
                    sequences.append(line.strip())
        f.close()

    # Create columns for data frame 
    fasta["Header"] = headers
    fasta["Isolate_Id"] = isolate_ids
    fasta["Isolate_Name"] = isolate_names
    fasta["Subtype"] = subtypes
    fasta["Segment"] = segments
    fasta["Collection_Date"] = collection_dates
    fasta["Sequence"] = sequences
    
    return fasta

# Function to get each unique animal listed so we can sort them

def sort_animals(fasta):
    isolate_names = fasta["Isolate_Name"]
    animal_list = []
    for name in isolate_names.values:
        animal = name.split("/")[1]
        animal_list.append(animal)

    unique_animals = list(set(animal_list))

    # Save the animals to a file so we can sort them
    return unique_animals

    # animals_df.to_csv(file_name) # We will have to sort manually, unfortunately :(

# Separate B3.13 and D1.1 in both the Excel file and the FASTA file

def separate_b313_d11_xls(metadata_raw, fasta):
    # Excel
    b313_xls = metadata_raw[metadata_raw["Genotype"] == "B3.13"]
    d11_xls = metadata_raw[metadata_raw["Genotype"] == "D1.1"]

    # print(d11_xls)

    # FASTA

    b313_mask = fasta['Isolate_Id'].isin(b313_xls['Isolate_Id'])
    d11_mask = fasta['Isolate_Id'].isin(d11_xls['Isolate_Id'])

    b313_fasta = fasta[b313_mask]
    d11_fasta = fasta[d11_mask]
    return b313_fasta, d11_fasta

# Fix animals in host type

def fix_animals(fasta, animals_ref):

    # If the animal is in a specific column of animals_ref, label host type as column name
    animal_list = [] # Find animals first
    for name in fasta["Isolate_Name"].values:
        animal = name.split("/")[1]
        animal_list.append(animal)

    animal_types = []
    for animal in animal_list: # Label each animal as a type
        # print(animal)
        if animal in animals_ref["avian"].values:
            animal_types.append("avian")
        elif animal in animals_ref["cattle"].values:
            animal_types.append("cattle")
        elif animal in animals_ref["feline"].values:
            animal_types.append("feline")
        elif animal in animals_ref["other_mammal"].values:
            animal_types.append("other_mammal")
        else: # If other
            animal_types.append("other")

    fasta["Host_Type"] = animal_types

# Separate FASTA files into 8 different files based on segment

def separate_fasta_by_seg(metadata, fasta, animals_df): #, b313_fasta, d11_fasta):

    # Dummy host type -- we'll actually add this in later
    fix_animals(fasta, animals_df)
    # b313_fasta["Host_Type"] = "other"
    # d11_fasta["Host_Type"] = "other"

    unique_segments = list(set(fasta["Segment"]))
    genotypes = ["B3.13", "D1.1"]
    # genotype_fastas = {"B3.13": b313_fasta, "D1.1": d11_fasta}
    # print(unique_segments)

    # “>Isolate_name|subtype|collection_date|host_type|genotype”

    segment_fastas = []
    for fasta_gen in genotypes: # .keys():
        for seg in unique_segments:

            xls = metadata[metadata["Genotype"] == fasta_gen]

            # print(d11_xls)

            # FASTA

            mask = fasta['Isolate_Id'].isin(xls['Isolate_Id'])

            fasta_seg_pre = fasta[mask]

            # fasta_seg = genotype_fastas[fasta_gen][genotype_fastas[fasta_gen]["Segment"] == seg]
            fasta_seg = fasta_seg_pre[fasta_seg_pre["Segment"] == seg]

            fasta_seg["Genotype"] = fasta_gen

            # Rename sequences 
            new_name = ">" + fasta_seg["Isolate_Name"] + "|" + fasta_seg["Subtype"] + "|" + fasta_seg["Collection_Date"] + "|" + fasta_seg["Host_Type"] + "|" + fasta_gen + "\n"
            fasta_seg["New_Name"] = new_name
            # print(fasta_seg["New_Name"])

            segment_fastas.append(fasta_seg)

    return segment_fastas, unique_segments
