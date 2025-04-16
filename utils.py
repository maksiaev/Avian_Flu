
import os
import glob 
import pandas as pd
import re
import requests
import xml.etree.ElementTree as ET
from collections import defaultdict 
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.common.action_chains import ActionChains
import time 

### GISAID functions ###

# Function to download data from GISAID
# Opening and downloading from the website

def open_gisaid(username, password, browser, sleep_time, start_date, end_date):

    if browser=="Firefox":
        # If you want to open Firefox
        driver = webdriver.Firefox()
    elif browser=="Chrome": # if Chrome...
        driver = webdriver.Chrome()
    else: # Edge, probably
        driver = webdriver.Edge()

    # How many seconds should pass between tasks
    sleep_time = int(sleep_time)

    # Requested URL
    driver.get("https://www.epicov.org/epi3/frontend#")

    # Wait for it to load, otherwise it won't work
    # driver.implicitly_wait(20)
    # username_field = wait.until(EC.element_to_be_clickable((By.NAME, 'login')))
    # password_field = wait.until(EC.element_to_be_clickable((By.NAME, 'password')))

    time.sleep(sleep_time)

    # Input username and password
    username_field = driver.find_element(By.NAME, "login")
    password_field = driver.find_element(By.NAME, "password")
    submit_button = driver.find_element(By.CLASS_NAME, "form_button_submit")
    username_field.send_keys(username)
    password_field.send_keys(password)


    # Wait for it to load
    time.sleep(sleep_time)
    # submit_button = wait.until(EC.element_to_be_clickable((By.CLASS_NAME, 'form_button_submit')))
    submit_button.click()

    # Wait for it to load again
    time.sleep(sleep_time)

    # Find the correct tab
    epiflu_link = driver.find_element(By.XPATH, "//*[contains(text(), 'EpiFlu™')]")
    # epiflu_link = wait.until(EC.visibility_of_element_located((By.XPATH, "//*[contains(text(), 'EpiFlu™')]")))
    epiflu_link.click()

    time.sleep(sleep_time)

    # Search tab 
    search_link = driver.find_element(By.XPATH, "//*[contains(text(), 'Search')]")
    search_link.click()

    time.sleep(sleep_time)

    # Type: A
    # Gives "internal server error" without action chains
    type_a = driver.find_element(By.XPATH, "//*[contains(@class, 'sys-form-filine-td')]//option[@value='A']")
    ActionChains(driver).move_to_element(type_a).pause(1).click(type_a).perform() 

    time.sleep(sleep_time)

    # H: 5
    h_5 = driver.find_element(By.XPATH, "//*[contains(@class, 'sys-form-filine-td')]//option[@value='5']") # The second multi-select
    ActionChains(driver).move_to_element(h_5).pause(1).click(h_5).perform() 

    time.sleep(sleep_time)

    # Find the ID for N

    # N: 1
    n_1 = driver.find_element(By.XPATH, "//*[contains(@class, 'sys-form-filine-td')][3]//option[@value='1']") # The third multi-select
    ActionChains(driver).move_to_element(n_1).pause(1).click(n_1).perform() 

    time.sleep(sleep_time)
    
    # Location: North America
    location_northa = driver.find_element(By.XPATH, "//*[contains(@class, 'sys-event-hook sys-fi-mark')]//option[@value='6440']")
    # location_northa = location_northa_1.location_once_scrolled_into_view
    driver.execute_script("arguments[0].scrollIntoView();", location_northa)
    ActionChains(driver).move_to_element(location_northa).pause(1).click(location_northa).perform()
    # location_northa.click()

    time.sleep(sleep_time)

    # Segments: PB2, PB1, PA, HA, NP, NA, MP, NS (all EXCEPT HE, P3)
    segment_list = ['PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'MP', 'NS']
    for segment in segment_list:
        xpath = "//*[contains(@class, 'sys-form-fi-cb sys-fi-mark')]//input[@value='" + segment + "']"
        pb2 = driver.find_element(By.XPATH, xpath)
        ActionChains(driver).move_to_element(pb2).pause(1).click(pb2).perform()

    # Start submission: 2023-03-18
    # End submission: 2025-03-31

    # Insert dates
    first_date = driver.find_element(By.XPATH, "//*[contains(text(),'Submission date from')]/ancestor::td/following-sibling::td//input[@class='sys-event-hook sys-fi-mark hasDatepicker']") # First date
    first_date.send_keys(start_date)

    last_date = driver.find_element(By.XPATH, "//*[contains(text(),'Submission date from')]/ancestor::td/following-sibling::td[3]//input[@class='sys-event-hook sys-fi-mark hasDatepicker']") # Second date
    last_date.send_keys(end_date)

    time.sleep(sleep_time)

    # Split search into batches by date if >= 10k sequences

    # Search

    search_button = driver.find_element(By.XPATH, "//*[contains(@class, 'buttons container-slot')]//button[@accesskey='g']")
    # search_button.click()
    ActionChains(driver).move_to_element(search_button).pause(1).click(search_button).perform() 

    # Choose all files -- unless the number of sequences is >10k, then do it in batches
    time.sleep(sleep_time)

    # # Double the sleep time to load
    # time.sleep(sleep_time)

    # Select all (at first)
    select_checkbox = driver.find_element(By.XPATH, "//*[contains(@class, 'yui-dt-first yui-dt-last')]//input[@type='checkbox']")
    select_checkbox.click()

    # Press the select button and select up to 10,000 entries. If there's more afterwards, download the first 10k and continue until we've reached the end.

    # Press the select button
    select_button = driver.find_element(By.XPATH, "//button[contains(text(), 'Select')]")
    ActionChains(driver).move_to_element(select_button).pause(1).click(select_button).perform()    

    time.sleep(sleep_time)

    # Get all the entries
    iframe = driver.find_element(By.NAME, "wjob")
    driver.switch_to.frame(iframe)
    entry_text = driver.find_element(By.XPATH, "//textarea")
    time.sleep(sleep_time)
    full_text = entry_text.text
    entries = re.split(r"[,]", full_text)

    # Get rid of empty strings
    for entry in entries:
        if entry == '':
            entries.remove(entry)

    # While there are >10k sequences, add each set of 10k to a list
    # 1 string ~= 8 sequences
    thresh = 500 # //10
    ten_thousands = []
    # entries = entries[thresh:]

    while len(entries) > thresh: 
        ten_k = entries[:thresh]
        ten_thousands.append(ten_k)
        entries = entries[thresh:]
        # print(len(entries))

    ten_thousands.append(entries) # if less than thresh to begin with, appends full batch. Else appends last batch. 

    # for t in ten_thousands:
    #     print(len(t))

    print(len(ten_thousands))

    # Go back and select the nth 10k from the entries frame

    # Clear the iframe
    entry_text.clear()

    for text in ten_thousands:
        for t in text:
            entry_text.send_keys(t)
        # Press OK
        ok_button = driver.find_element(By.XPATH, "//button[contains(text(), 'OK')]")
        ok_button.click()

        time.sleep(sleep_time)
        # Press OK again
        ok_button2 = driver.find_element(By.XPATH, "//*[contains(@class, 'yui-panel-container shadow')]//button[contains(text(), 'OK')]")
        ok_button2.click()

        time.sleep(sleep_time)

        #Press download
        driver.switch_to.default_content()
        download_button = driver.find_element(By.XPATH, "//*[contains(@class, 'buttons container-slot')]//button[contains(text(), 'Download')]")
        ActionChains(driver).move_to_element(download_button).pause(1).click(download_button).perform()    

        time.sleep(sleep_time) 

        # Download metadata as XLS
        iframe_download = driver.find_element(By.NAME, "downl")
        driver.switch_to.frame(iframe_download)
        time.sleep(sleep_time)
        download_button = driver.find_element(By.XPATH, "//*[contains(@class, 'sys-component-slot')]//button[contains(text(), 'Download')]")
        ActionChains(driver).move_to_element(download_button).pause(1).click(download_button).perform()   
        # download_button.click()

        time.sleep(sleep_time * 10) # Downloads can take a while

        #Press download
        driver.switch_to.default_content()
        download_button = driver.find_element(By.XPATH, "//*[contains(@class, 'buttons container-slot')]//button[contains(text(), 'Download')]")
        ActionChains(driver).move_to_element(download_button).pause(1).click(download_button).perform()    

        time.sleep(sleep_time)

        # Download segment sequences as DNA FASTA file
        iframe_download = driver.find_element(By.NAME, "downl")
        driver.switch_to.frame(iframe_download)
        time.sleep(sleep_time)
        dna = driver.find_element(By.XPATH, "//input[@value='dna']")
        ActionChains(driver).move_to_element(dna).pause(1).click(dna).perform()    

        time.sleep(sleep_time)

        # Rename
        name_text = driver.find_element(By.XPATH, "//input[@type='text']")
        name_text.clear()
        t = "Isolate ID | Isolate name | Type | Segment | Collection date"
        name_text.send_keys(t)
        time.sleep(sleep_time)
        # Segments: PB2, PB1, PA, HA, NP, NA, MP, NS (all EXCEPT HE, P3)
        segment_list = ['PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'MP', 'NS']
        for segment in segment_list:
            xpath = "//input[@value='" + segment + "']"
            pb2 = driver.find_element(By.XPATH, xpath)
            ActionChains(driver).move_to_element(pb2).pause(1).click(pb2).perform()
        # Assume both checkboxes for spaces and underscores are already checked -- if not, come back to this
        download_button = driver.find_element(By.XPATH, "//*[contains(@class, 'sys-component-slot')]//button[contains(text(), 'Download')]")
        ActionChains(driver).move_to_element(download_button).pause(1).click(download_button).perform()   

        time.sleep(sleep_time)

        go_back = driver.find_element(By.XPATH, "//button[contains(text(), 'Go back')]")
        ActionChains(driver).move_to_element(go_back).pause(1).click(go_back).perform()   

        time.sleep(sleep_time)

        driver.switch_to.default_content()

        time.sleep(sleep_time * 10)

        # Press the select button
        select_button = driver.find_element(By.XPATH, "//*[contains(@class, 'buttons container-slot')]//button[contains(text(), 'Select')]")
        ActionChains(driver).move_to_element(select_button).pause(1).click(select_button).perform()

        time.sleep(sleep_time)  

        # Clear the iframe
        iframe = driver.find_element(By.NAME, "wjob")
        driver.switch_to.frame(iframe)
        entry_text = driver.find_element(By.XPATH, "//textarea")
        entry_text.clear()

    # Close browser

    driver.close()

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
            # print(line)
            if line[0] == ">": # If it's a header
                if line[1:].strip() not in headers: # And the previous line is not a header we've seen before
                    header = line[1:].strip() # Remove the ">"
                    # print(header)
                    split_header = header.split("|")
                    # print(split_header)
                    headers.append(header) 
                    isolate_ids.append(split_header[0])
                    isolate_names.append(split_header[1]) # We'll need to extract data from this too
                    # print(split_header[2].split("_")[-1])
                    subtypes.append(split_header[2].split("_")[-1])  # Get only H5N1
                    segments.append(split_header[3])
                    if split_header[4] == "2024-01-01":
                        collection_dates.append("2024") # No samples were collected 1/1/2024, these are all unknown 
                    elif split_header[4] == "2025-01-01":
                        collection_dates.append("2025")
                    else: 
                        collection_dates.append(split_header[4])
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

# Separate B3.13 and D1.1 in both the Excel file and the FASTA file

# def separate_b313_d11_xls(metadata_raw, fasta):
#     # Excel
#     b313_xls = metadata_raw[metadata_raw["Genotype"] == "B3.13"]
#     d11_xls = metadata_raw[metadata_raw["Genotype"] == "D1.1"]

#     # print(d11_xls)

#     # FASTA

#     b313_mask = fasta['Isolate_Id'].isin(b313_xls['Isolate_Id'])
#     d11_mask = fasta['Isolate_Id'].isin(d11_xls['Isolate_Id'])

#     b313_fasta = fasta[b313_mask]
#     d11_fasta = fasta[d11_mask]
#     return b313_fasta, d11_fasta

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
        if animal in animals_ref["avian"].values:
            animal_types.append("avian")
        elif animal in animals_ref["cattle"].values:
            animal_types.append("cattle")
        elif animal in animals_ref["feline"].values:
            animal_types.append("feline")
        elif animal in animals_ref["other_mammal"].values:
            animal_types.append("other_mammal")
        elif animal in animals_ref["human"].values:
            animal_types.append("human")
        else: # If other
            animal_types.append("other")

    fasta["Host_Type"] = animal_types

    return fasta

# Separate FASTA files into 8 different files based on segment

def separate_fasta_by_seg(metadata, fasta, animals_df, genotypes): #, b313_fasta, d11_fasta):

    fasta = fix_animals(fasta, animals_df)
    # Dummy host type -- we'll actually add this in later
    # b313_fasta["Host_Type"] = "other"
    # d11_fasta["Host_Type"] = "other"

    unique_segments = list(set(fasta["Segment"]))
    # genotypes = ["B3.13", "D1.1"]
    # genotype_fastas = {"B3.13": b313_fasta, "D1.1": d11_fasta}

    # “>Isolate_name|subtype|collection_date|host_type|genotype”

    segment_fastas = []
    for fasta_gen in genotypes: # .keys():
        for seg in unique_segments:

            xls = metadata[metadata["Genotype"].apply(lambda x: x.split(" ")[0]) == fasta_gen]

            # print("XLS: ", metadata["Genotype"])

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



### Andersen Lab functions ###

# Rename host type

def sort_animals_andersen(metadata):
    host_names = metadata["Host"]
    animal_list = []
    for name in host_names.values:
        animal_low = name.lower()
        # animal = animal_low.replace(" ", "_") # replace spaces with underscores
        animal_list.append(animal_low)

    unique_animals = list(set(animal_list))

    # Save the animals to a file so we can sort them
    return unique_animals

# Fix animals

def fix_animals_andersen(metadata, animals_ref):

    # If the animal is in a specific column of animals_ref, label host type as column name
    animal_list = [] # Find animals first
    for name in metadata["Host"].values:
        name = name.lower()
        # name = name.replace(" ", "_")
        animal_list.append(name)

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

    metadata["Host_Type"] = animal_types

# Get collection date from GenBank eutils 

def search_collection_date(biosample, metadata_genbank):

    print(biosample)

    try:

        # Avoid spamming the server
        time.sleep(2)
    
        base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
        search_url = base_url + "esearch.fcgi?db=biosample&term=" + biosample +"&usehistory=y&api_key=2cbaf77ac9ec5ae7844ea350076ae6d56809"

        # Get Biosample ID from search_url
        output = requests.get(search_url)
        xml = output.content
        root = ET.fromstring(xml)
        sample_id = root.find("./IdList/Id").text

        biosample_url = base_url + "elink.fcgi?dbfrom=biosample&db=nuccore&id=" + sample_id + "&cmd=neighbor_history&api_key=2cbaf77ac9ec5ae7844ea350076ae6d56809"
        
        # Get Nucleotide ID from biosample_url
        output = requests.get(biosample_url)
        xml = output.content
        root = ET.fromstring(xml)
        query_key = root.find(".//QueryKey").text
        web_env = root.find(".//WebEnv").text

        nucleotide_url = base_url + "esummary.fcgi?db=nuccore&query_key=" + query_key + "&WebEnv=" + web_env + "&version=2.0&api_key=2cbaf77ac9ec5ae7844ea350076ae6d56809"

        output = requests.get(nucleotide_url) 
        xml = output.content
        root = ET.fromstring(xml)

        # Grab collection date at the end of the sub name
        collection_date = root.find(".//SubName").text.split("|")[-1]

        return collection_date
    
    except:
        print("Unable to find collection date.")

        if len(metadata_genbank[metadata_genbank["BioSample"] == biosample]["Collection_Date"]) > 0: # If a year exists
            collection_date = metadata_genbank[metadata_genbank["BioSample"] == biosample]["Collection_Date"].values[0]
        else:
            collection_date = float('nan') 

        return collection_date
    
def create_dataframes(directory):
    dfs_gisaid = defaultdict(list)
    for dirpath, dirs, files in os.walk(directory): # Find the fasta file
        for file in files:
            file_name = os.path.join(dirpath, file) # Get file name
            # print(file_name)
            gisaid_df = pd.DataFrame()
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
                        other = ""
                        for d in digits:
                            # print(d)
                            if len(d) == 6 and d.isnumeric(): # If it's just digits and not one of those weird isolates
                                isolate = d + "-"
                            elif len(d) == 3 and d.isnumeric():
                                isolate = isolate + d
                            elif d.isnumeric() == False: # If it's a weird isolate
                                other = d + "-"
                            else: 
                                other = other + d
                        # Now add to list to check in Andersen files without doing wild for loops
                        if len(isolate) == 10: # If this is a correctly formatted isolate
                            # isolates.append(isolate)
                            # All headers are followed by sequences
                            isolate_partial.append(isolate)
                        else: # If this is some other isolate
                            isolate_partial.append(other)
                    elif line == "nan\n":
                        print(directory) # Some headers in the Andersen files don't exist 
                        isolate_partial.append(float('nan'))
                        full_header.append(line) # Sorry :/
                    else: # It's a sequence
                        sequence.append(line)
                    
                gisaid_df["isolate_partial"] = isolate_partial
                gisaid_df["full_header"] = full_header
                # print(len(sequence))
                gisaid_df["sequence"] = sequence
                dfs_gisaid[file_name.split("/")[-1][:-6]].append(gisaid_df)
    return dfs_gisaid