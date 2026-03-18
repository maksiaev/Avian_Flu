# Housekeeping
import os
import pandas as pd
import re
import dateutil 
from collections import defaultdict 
import time 

# ## NCBI Virus functions ###

# ## Andersen Lab functions ###

# Rename host type
def sort_animals_andersen(metadata):
    host_names = metadata["Host"]
    animal_list = []
    for name in host_names.values:
        try:
            if name.replace("-", " ").replace(" ", "").replace(".", "").replace(",", "").isalpha(): # If all characters are alphanumeric
                animal_low = name.lower() # make the name lowercase
                animal_list.append(animal_low)
        except:
            animal_list.append("unknown")

    unique_animals = list(set(animal_list))

    # Save the animals to a file so we can sort them (outside of this function)
    return unique_animals

# Fix animals
def fix_animals_andersen(metadata, animals_ref):

    # If the animal is in a specific column of animals_ref, label host type as column name
    animal_list = [] # Find animals first
    for name in metadata["Host"].values:
        try:
            name = name.lower()
            animal_list.append(name)
        except:
            animal_list.append("unknown")

    animal_types = []
    for animal in animal_list: # Label each animal as a type
        # print(animal)
        if animal in animals_ref["wild_avian"].values:
            animal_types.append("wild_avian")
        elif animal in animals_ref["domestic_avian"].values:
            animal_types.append("domestic_avian")
        elif animal in animals_ref["cattle"].values:
            animal_types.append("cattle")
        elif animal in animals_ref["feline"].values:
            animal_types.append("feline")
        elif animal in animals_ref["human"].values:
            animal_types.append("human")
        elif animal in animals_ref["other_mammal"].values:
            animal_types.append("other_mammal")
        elif animal in animals_ref["pet_food"].values:
            animal_types.append("pet_food")
        else: # If other
            animal_types.append("other")

    metadata["Host_Type"] = animal_types
    # metadata["Host"] = metadata["Host"].apply(lambda x: x.replace(" ", "_")) # Avoid space issues

def create_dataframes(directory):
    dfs = defaultdict(list)
    for dirpath, dirs, files in os.walk(directory): # Find the fasta file
        for file in files:
            file_name = os.path.join(dirpath, file) # Get file name
            # print(file_name)
            df = pd.DataFrame()
            with open(file_name) as f:
                lines = f.readlines()
                isolate_partial = []
                full_header = []
                sequence = []
                # Some lines start with 25_, others 25-. This shouldn't matter, but split on "_" first
                for num, line in enumerate(lines):
                    if line[0] == ">": # If it's a header
                        full_header.append(line)
                        full = line.split("/")[3] # Get the isolate
                        partial = full.split("_")[-1] # If 25_, get the last bit
                        digits = partial.split("-")
                        isolate = ""
                        # other = ""
                        for d in digits:
                            # print(d)
                            if len(d) == 6 and d.isnumeric(): # If it's just digits and not one of those weird isolates
                                isolate = d + "-"
                            elif len(d) == 3 and d.isnumeric():
                                isolate = isolate + d
                            # elif d.isnumeric() == False: # If it's a weird isolate
                            #     other = d + "-"
                            # else: # If it's a weird isolate
                            #     other = other + d
                        # Now add to list to check in Andersen files without doing wild for loops
                        if len(isolate) == 10: # If this is a correctly formatted isolate
                            # isolates.append(isolate)
                            # All headers are followed by sequences
                            isolate_partial.append(isolate)
                        else: # If this is some other isolate
                            isolate_partial.append(full)
                    elif line == "nan\n":
                        print(directory) # Some headers in the Andersen files don't exist 
                        isolate_partial.append(float('nan'))
                        full_header.append(line) # Sorry :/
                    else: # It's a sequence
                        sequence.append(line)
                    
                df["isolate_partial"] = isolate_partial
                df["full_header"] = full_header
                # print(len(sequence))
                df["sequence"] = sequence
                key_name = (file_name.split("/")[-1][:-6]).split("_")[0] + "_" + (file_name.split("/")[-1][:-6]).split("_")[1]
                print(key_name)
                dfs[key_name].append(df)
            f.close()
        break # Don't go into subdirectories 
    return dfs

# Function to prepare dataframes from complete fastas
def fasta_df_complete(file_name, state_ref):

    fasta = pd.DataFrame()
    headers = []
    isolate_ids = []
    isolate_partials = []
    isolate_names = []
    subtypes = []
    # segments = []
    collection_dates = []
    sequences = []
    host_types = []
    species = []
    identifiers = []
    genotypes = []
    with open(file_name) as f:
        lines = f.readlines()
        for num, line in enumerate(lines):
            # print(line)
            if line[0] == ">": # If it's a header
                if line[1:].strip() not in headers: # And the previous line is not a header we've seen before
                    header = line[1:].strip() # Remove the ">"
                    # print(header)
                    try:
                        split_header = header.split("|")
                        if len(header.split("|")) > 6 and len(header.split("|")) < 8:
                            identifier = header.split("|")[0]
                            identifiers.append(identifier)
                            split_first_header = split_header[1].split("/")
                        elif len(header.split("|")) >= 8:
                            identifier = header.split("|")[0]
                            identifiers.append(identifier)
                            split_first_header = split_header[1].split("/")
                        else:
                            identifier = "unknown"
                            identifiers.append(identifier)
                            split_first_header = split_header[0].split("/")
                        # print(split_first_header)
                        # print(split_header)
                        headers.append(header) 
                        isolate_ids.append(split_first_header[3])
                        isolate_partials.append(partial_isolate(split_first_header[3]))
                        isolate_names.append("/".join(split_first_header)) # We'll need to extract data from this too
                        # print(split_header[2].split("_")[-1])
                        subtypes.append(split_header[-5].split("_")[-1])  # Get only H5N1
                        genotypes.append(split_header[-1])
                        # segments.append(split_header[-2]) # .split("|")[-1])
                        host_types.append(split_header[-2])
                        species.append(split_first_header[1])
                        # if split_header[4] == "2024-01-01":
                        #     collection_dates.append("2024") # No samples were collected 1/1/2024, these are all unknown 
                        # elif split_header[4] == "2025-01-01":
                        #     collection_dates.append("2025")
                        # else: 
                        collection_dates.append(split_header[-3])
                        # collection_dates.append(split_header[-1].split("_")[-1])
                        if num < len(lines): # If we're not at the last line
                            # for i, l in enumerate(lines[num + 1:]):
                            i = num
                            sequence = ""
                            # print(lines[i])
                            # print(lines[i + 1])
                            while i < len(lines) - 1 and lines[i + 1][0] != ">": # While the next line is part of a sequence
                                sequence = sequence + lines[i + 1].strip()
                                i += 1
                            sequences.append(sequence) # Add next line to sequences
                    except:
                        if header in headers:
                            headers.remove(header)
                        if identifier in identifiers:
                            identifiers.remove(identifier)
                        print(header)
                        continue
        f.close()

    # Create columns for data frame 
    fasta["Header"] = headers
    fasta["Isolate_Id"] = isolate_ids
    fasta["Isolate_Name"] = isolate_names
    fasta["Subtype"] = subtypes
    fasta["Partials"] = isolate_partials
    # fasta["Segment"] = segments
    # Geo_Location is more complicated
    fasta["Location"] = fasta["Isolate_Name"].apply(lambda x: x.split("/")[-3])
    fasta["Geo_Location"] = fasta["Location"].apply( # lambda x: state_ref.loc[state_ref["Abbreviation"] == x.split("/")[2], 'Country'].iloc[0] + "-" + x.split("/")[2] if x.split("/")[2] in state_ref["Abbreviation"].values else state_ref.loc[state_ref["State"] == x.split("/")[2].replace("_", " "), 'Country'].iloc[0] + "-" + state_ref.loc[state_ref["State"] == x.split("/")[2].replace("_", " "), 'Abbreviation'].iloc[0] if x.split("/")[2].replace("_", " ") in state_ref["State"].values else x.split("/")[2].replace(": ", "-"))

                                                    lambda x: 
                                                    # If "x" has the state abbreviation (e.g. "MD")
                                                    state_ref.loc[state_ref["Abbreviation"].str.contains('|'.join(x.replace(": ", ",").replace(" ", "_").split(',')), regex=True), 'Country'].iloc[0] 
                                                    + "-" + 
                                                    x.split(" ")[-1]
                                                    if state_ref["Abbreviation"].str.contains("|".join((x.replace(": ", ",").replace(" ", "_").split(','))), regex=True).any()
                                                    # If "x" has the full state name (e.g. "Maryland")
                                                    else state_ref.loc[state_ref['State'].str.contains('|'.join(x.replace(": ", ",").replace(" ", "_").split(',')), regex=True), 'Country'].iloc[0]
                                                    + "-" + 
                                                    state_ref.loc[state_ref['State'].str.contains('|'.join(x.replace(": ", ",").replace(" ", "_").split(',')), regex=True), 'Abbreviation'].iloc[0] 
                                                    if state_ref["State"].str.contains("|".join((x.replace(": ", ",").replace(" ", "_").split(','))), regex=True).any()
                                                    # If "x" has neither the state abbreviation nor the full state name nor is "USA"
                                                    else 
                                                    x
                                                    )
    fasta["Date Collected"] = collection_dates
    fasta["Species"] = species
    fasta["Host_Type"] = host_types
    fasta["Genotype"] = genotypes
    fasta["Sequence"] = sequences
    # if len(identifiers) == len(fasta):
    fasta["Identifier"] = identifiers
    
    return fasta

def df_from_fasta(file_name):
    df = pd.DataFrame()
    with open(file_name) as f:
        lines = f.readlines()
        full_header = []
        sequences = []
        sequence = ""
        for num, line in enumerate(lines):
            if line[0] == ">" and num < 1: # If it's the first header
                full_header.append(line)
            elif line[0] == ">" and num >= 1: # If it's a subsequent header
                full_header.append(line)
                sequences.append(sequence)
                sequence = "" # Erase sequence 
            else: # Else it's a sequence
                sequence += line
            if num == len(lines) - 1: # If it's the last part of the last sequence
                sequences.append(sequence)
            
        df["full_header"] = full_header
        df["sequence"] = sequences
    f.close()
    return df

def df_to_fasta(fasta, file_name, output_path):

    output_file = open(output_path + file_name, "a")

    for index, row in fasta.iterrows():
        name = fasta.loc[index, "full_header"]
        sequence = fasta.loc[index, "sequence"]
    # First is header, second is sequence
        output_file.write(name + "\n")
        output_file.write(sequence + "\n")
    output_file.close()

def relabel_animals(fasta, animals_ref):

    names = []
    for name in fasta["full_header"].values:
        animal = name.split("/")[1]
        animal = animal.lower()
        host_type = name.split("|")[-2]

        if animal in animals_ref["wild_avian"].values:
            new_name = name.replace(host_type, "wild_avian")
            names.append(new_name)
        elif animal in animals_ref["domestic_avian"].values:
            new_name = name.replace(host_type, "domestic_avian")
            names.append(new_name)
        elif animal in animals_ref["cattle"].values:
            new_name = name.replace(host_type, "cattle")
            names.append(new_name)
        elif animal in animals_ref["feline"].values:
            new_name = name.replace(host_type, "feline")
            names.append(new_name)
        elif animal in animals_ref["other_mammal"].values:
            new_name = name.replace(host_type, "other_mammal")
            names.append(new_name)
        elif animal in animals_ref["human"].values:
            new_name = name.replace(host_type, "human")
            names.append(new_name)
        elif animal in animals_ref["pet_food"].values:
            new_name = name.replace(host_type, "pet_food")
            names.append(new_name)
        else: # If other
            new_name = name.replace(host_type, "other")
            names.append(new_name)

    fasta["full_header"] = names

    return fasta

# ## GISAID functions ###

# Function to get each unique animal listed so we can sort them

def sort_animals(fasta):
    isolate_names = fasta["Isolate_Name"]
    animal_list = []
    for name in isolate_names.values:
        # print(name)
        try:
            animal = name.split("/")[1]
            animal_low = animal.lower()
            animal_list.append(animal_low)
        except:
            continue

    unique_animals = list(set(animal_list))

    # Save the animals to a file so we can sort them
    return unique_animals

    # animals_df.to_csv(file_name) # We will have to sort manually, unfortunately :(

# Fix animals in host type

def fix_animals(fasta, animals_ref):

    # If the animal is in a specific column of animals_ref, label host type as column name
    animal_list = [] # Find animals first
    for name in fasta["Isolate_Name"].values:
        try:
            animal = name.split("/")[1]
            animal_low = animal.lower()
        except:
            animal_low = "unknown"
        animal_list.append(animal_low)


    animal_types = []
    for animal in animal_list: # Label each animal as a type
        # print(animal)
        if animal in animals_ref["wild_avian"].values:
            animal_types.append("wild_avian")
        elif animal in animals_ref["domestic_avian"].values:
            animal_types.append("domestic_avian")
        elif animal in animals_ref["cattle"].values:
            animal_types.append("cattle")
        elif animal in animals_ref["feline"].values:
            animal_types.append("feline")
        elif animal in animals_ref["other_mammal"].values:
            animal_types.append("other_mammal")
        elif animal in animals_ref["human"].values:
            animal_types.append("human")
        elif animal in animals_ref["pet_food"].values:
            animal_types.append("pet_food")
        else: # If other
            animal_types.append("other")

    fasta["Host_Type"] = animal_types

    return fasta

def partial_isolate(id):

    partial = id.replace("_", "-") # [-1] # If 25_, get the last bit
    partial = partial.replace("-original", "") # If -original suffix, remove
    partial = partial.replace("G", "-0") # If 25G, fix
    digits = partial.split("-")
    # Build partial isolates
    isolate = ""
    # other = ""
    for d in digits:
        # print(d)
        if len(d) == 2 and d.isnumeric(): # The 25 bit
            # isolate = d + "-"
            pass
        if len(d) == 6 and d.isnumeric(): # If it's just digits and not one of those weird isolates
            isolate = d + "-"
        if len(d) == 3 and d.isnumeric():
            isolate = isolate + d
        # elif d.isnumeric() == False: # If it's a weird isolate
        #     other = d + "-"
        # else: 
        #     other = partial 
    # Now add to list to check in Andersen files without doing wild for loops
    if len(isolate) == 10: # If this is a correctly formatted isolate
        # isolates.append(isolate)
        # All headers are followed by sequences
        partial_isolate = isolate
    else: # If this is some other isolate
        partial_isolate = partial

    return partial_isolate

def metadata_from_fasta(file_name):
    fasta_df = df_from_fasta(file_name)
    try:
        fasta_df["Accession"] = fasta_df["full_header"].apply(lambda x: x.split("|")[0].replace(">", ""))
        fasta_df["Host"] = fasta_df["full_header"].apply(lambda x: "human" if "human" in x else x.split("/")[1] if "/" in x else "unknown") #  if len(x.split("/")) > 1 else "unknown")
        fasta_df["Isolate"] = fasta_df["full_header"].apply(lambda x: x.split("/")[2] if "human" in x else x.split("/")[3] if "/" in x else "unknown") #  if len(x.split("/")) > 2 else "unknown")
        fasta_df["Subtype"] = fasta_df["full_header"].apply(lambda x: x.split("|")[-5])
        fasta_df["Geo_Location"] = fasta_df["full_header"].apply(lambda x: x.split("|")[-4])
        fasta_df["Collection_Date"] = fasta_df["full_header"].apply(lambda x: x.split("|")[-3])
        fasta_df["Host_Type"] = fasta_df["full_header"].apply(lambda x: x.split("|")[-2])
        fasta_df["Genotype"] = fasta_df["full_header"].apply(lambda x: x.split("|")[-1])
    except:
        print(fasta_df)
    return fasta_df

# Function to get metadata
def separate_fasta_by_segs(metadata, fasta, animals_df, genotypes): #, b313_fasta, d11_fasta):

    fasta = fix_animals(fasta, animals_df) # Fix animals first

    unique_segments = list(set(fasta["Segment"])) # Get list of segments

    # “>EPI_ID|Isolate_name|subtype|collection_date|host_type|genotype”

    segment_fastas = [] # Get a list of fastas, separated by segment
    for genotype in genotypes: # .keys(): # For each genotype
        print(genotype)
        for seg in unique_segments: # For each segment

            xls = metadata # [metadata["Genotype"].apply(lambda x: x.split(" ")[0]) == genotype] # Get only the metadata corresponding to that genotype
            xls = xls.rename(columns={"Isolate_Id":"Identifier"})
            print(xls)
            # print("XLS: ", metadata["Genotype"])
            # print(xls["Clade"])

            # fasta_seg_pre = fasta[fasta["Identifier"].isin(xls['Isolate_Id'])] # Get only the identifiers (Isolate_Id) that are left after metadata is filtered for genotype
            fasta_seg_pre = fasta.merge(xls, how="left", on="Identifier")
            print(fasta_seg_pre.columns)
            # break 

            # fasta_seg = genotype_fastas[fasta_gen][genotype_fastas[fasta_gen]["Segment"] == seg]
            fasta_seg = fasta_seg_pre[fasta_seg_pre["Segment"] == seg]

            fasta_seg["Genotype"] = genotype

            # Rename sequences 
            new_name = ">" + fasta_seg["Identifier"] + "|" + fasta_seg["Isolate_Name_x"] + "|" + fasta_seg["Subtype_y"] + "|" + fasta_seg["Location"].apply(lambda x: x.split(" / ")[-1]) + "|" + fasta_seg["Date Collected"].apply(lambda x: str(dateutil.parser.parse(x, default=dateutil.parser.parse("2000-01-01")).year) if dateutil.parser.parse(x, default=dateutil.parser.parse("2000-01-01")).month == dateutil.parser.parse("2000-01-01").month and dateutil.parser.parse(x, default=dateutil.parser.parse("2000-01-01")).day == dateutil.parser.parse("2000-01-01").day else x) + fasta_seg["Host_Type"] #.apply(lambda x: "" if x != "human" else "|human")
            fasta_seg["full_header"] = new_name
            # print(fasta_seg["New_Name"])

            segment_fastas.append(fasta_seg)
            print(fasta_seg)

            segment_fastas.append(fasta_seg)

    if "Unassigned" in genotypes: # if we care about the unassigned sequences
        # Unassigned genotypes
        for seg in unique_segments:
            xls = metadata[metadata["Genotype"].str.contains("Not")] # Get only the metadata corresponding to that genotype

            xls = xls.rename(columns={"Isolate_Id":"Identifier"})
            print(xls)
            # print("XLS: ", metadata["Genotype"])
            # print(xls["Clade"])

            # fasta_seg_pre = fasta[fasta["Identifier"].isin(xls['Isolate_Id'])] # Get only the identifiers (Isolate_Id) that are left after metadata is filtered for genotype
            fasta_seg_pre = fasta.merge(xls, how="left", on="Identifier")
            print(fasta_seg_pre.columns)
            # break 

            # fasta_seg = genotype_fastas[fasta_gen][genotype_fastas[fasta_gen]["Segment"] == seg]
            fasta_seg = fasta_seg_pre[fasta_seg_pre["Segment"] == seg]

            fasta_seg["Genotype"] = "Unassigned"

            # Rename sequences 
            new_name = ">" + fasta_seg["Identifier"] + "|" + fasta_seg["Isolate_Name_x"] + "|" + fasta_seg["Subtype_y"] + "|" + fasta_seg["Location"].apply(lambda x: x.split(" / ")[-1]) + "|" + fasta_seg["Date Collected"].apply(lambda x: str(dateutil.parser.parse(x, default=dateutil.parser.parse("2000-01-01")).year) if dateutil.parser.parse(x, default=dateutil.parser.parse("2000-01-01")).month == dateutil.parser.parse("2000-01-01").month and dateutil.parser.parse(x, default=dateutil.parser.parse("2000-01-01")).day == dateutil.parser.parse("2000-01-01").day else x) + fasta_seg["Host_Type"] #.apply(lambda x: "" if x != "human" else "|human")
            fasta_seg["full_header"] = new_name
            # print(fasta_seg["New_Name"])

            segment_fastas.append(fasta_seg)
            print(fasta_seg)

            segment_fastas.append(fasta_seg)

    return segment_fastas, unique_segments

def fasta_df(file_name, state_ref):
    fasta = pd.DataFrame()
    headers = []
    isolate_ids = []
    isolate_names = []
    subtypes = []
    segments = []
    collection_dates = []
    sequences = []
    # host_types = []
    species = []
    identifiers = []
    locations = []
    # genotypes = []
    with open(file_name) as f:
        lines = f.readlines()
        for num, line in enumerate(lines):
            # print(line)
            if line[0] == ">": # If it's a header
                if line[1:].strip() not in headers: # And the previous line is not a header we've seen before
                    header = line[1:].strip() # Remove the ">"
                    # print(header)
                    try:
                        split_header = header.split("|")
                        if len(header.split("|")) > 4:
                            identifier = header.split("|")[0]
                            identifiers.append(identifier)
                            split_first_header = split_header[1].split("/")
                        # else:
                        #     identifiers.append("unknown")
                        #     split_first_header = split_header[0].split("/")
                        # print(split_first_header)
                        # print(split_header)
                        headers.append(header) 
                        isolate_ids.append(split_first_header[-2]) # if "Catalonia" not in split_first_header[1] and "Navarra" not in split_first_header[1] else split_first_header[2])
                        locations.append(split_first_header[-3])
                        isolate_names.append(split_header[-4]) # We'll need to extract data from this too
                        # print(split_header[2].split("_")[-1])
                        subtypes.append(split_header[-3].split("_")[-1])  # Get only H5N1
                        # genotypes.append(split_header[-1])
                        segments.append(split_header[-2]) # .split("|")[-1])
                        # host_types.append(split_header[-2])
                        species.append(split_first_header[1])
                        # if split_header[4] == "2024-01-01":
                        #     collection_dates.append("2024") # No samples were collected 1/1/2024, these are all unknown 
                        # elif split_header[4] == "2025-01-01":
                        #     collection_dates.append("2025")
                        # else: 
                        collection_dates.append(split_header[-1])
                        # collection_dates.append(split_header[-1].split("_")[-1])
                        if num < len(lines): # If we're not at the last line
                            # for i, l in enumerate(lines[num + 1:]):
                            i = num
                            sequence = ""
                            # print(lines[i])
                            # print(lines[i + 1])
                            while i < len(lines) - 1 and lines[i + 1][0] != ">": # While the next line is part of a sequence
                                sequence = sequence + lines[i + 1].strip()
                                i += 1
                            sequences.append(sequence) # Add next line to sequences
                    except:
                        headers.remove(header)
                        identifiers.remove(identifier)
                        print(header)
                        continue
        f.close()

    
    # Create columns for data frame 
    fasta["Header"] = headers
    fasta["Isolate_Id"] = isolate_ids
    fasta["Isolate_Name"] = isolate_names
    fasta["Subtype"] = subtypes
    fasta["Segment"] = segments
    fasta["Location_Header"] = locations
    # Geo_Location is more complicated
    try:
        fasta["Geo_Location"] = fasta["Location_Header"].apply( # lambda x: state_ref.loc[state_ref["Abbreviation"] == x.split("/")[2], 'Country'].iloc[0] + "-" + x.split("/")[2] if x.split("/")[2] in state_ref["Abbreviation"].values else state_ref.loc[state_ref["State"] == x.split("/")[2].replace("_", " "), 'Country'].iloc[0] + "-" + state_ref.loc[state_ref["State"] == x.split("/")[2].replace("_", " "), 'Abbreviation'].iloc[0] if x.split("/")[2].replace("_", " ") in state_ref["State"].values else x.split("/")[2].replace(": ", "-"))

                                                    lambda x: 
                                                    # If "x" has the state abbreviation (e.g. "MD")
                                                    state_ref.loc[state_ref["Abbreviation"].str.contains('|'.join(x.replace(": ", ",").replace(" ", "_").split(',')), regex=True), 'Country'].iloc[0] 
                                                    + "-" + 
                                                    x.split(" ")[-1]
                                                    if state_ref["Abbreviation"].str.contains("|".join((x.replace(": ", ",").replace(" ", "_").split(','))), regex=True).any()
                                                    # If "x" has the full state name (e.g. "Maryland")
                                                    else state_ref.loc[state_ref['State'].str.contains('|'.join(x.replace(": ", ",").replace(" ", "_").split(',')), regex=True), 'Country'].iloc[0]
                                                    + "-" + 
                                                    state_ref.loc[state_ref['State'].str.contains('|'.join(x.replace(": ", ",").replace(" ", "_").split(',')), regex=True), 'Abbreviation'].iloc[0] 
                                                    if state_ref["State"].str.contains("|".join((x.replace(": ", ",").replace(" ", "_").split(','))), regex=True).any()
                                                    # If "x" has neither the state abbreviation nor the full state name nor is "USA"
                                                    else 
                                                    x
                                                    )
        fasta["Geo_Location"] = fasta["Geo_Location"].apply(lambda x: x.split("-")[0] if x.split("-")[-1] == "" or x.split("-")[-1] == x.split("-")[0] else x)
    except:
        print("Geo Location not found.")
        fasta["Geo_Location"] = fasta["Location_Header"]
    fasta["Date Collected"] = collection_dates
    fasta["Date Collected"] = fasta["Date Collected"].apply(lambda x: str(dateutil.parser.parse(x, default=dateutil.parser.parse("2000-01-01")).year) if dateutil.parser.parse(x, default=dateutil.parser.parse("2000-01-01")).month == dateutil.parser.parse("2000-01-01").month and dateutil.parser.parse(x, default=dateutil.parser.parse("2000-01-01")).day == dateutil.parser.parse("2000-01-01").day else x)
    fasta["Species"] = species
    # fasta["Host_Type"] = host_types
    # fasta["Genotype"] = genotypes
    fasta["Sequence"] = sequences
    # if len(identifiers) == len(fasta):
    fasta["Identifier"] = identifiers
        
    return fasta
