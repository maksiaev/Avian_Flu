from Bio import Phylo
from Bio import Nexus
from Bio import SeqIO
import os
import pandas as pd
import numpy as np

os.chdir("C:/Users/maksiaevai.NCBI_NT/Documents/Avian_Flu_Files/Trees/D1_3/")

# records = SeqIO.parse("D1_3_11-01-2021--07-04-2025_concat_318_tree_nexus.nex", "nexus")
# count = SeqIO.write(records, "D1_3_11-01-2021--07-04-2025_concat_318_tree_nexus_relabeled.nexus", "nexus")

# Load the Nexus file
# Replace 'your_nexus_file.nex' with the actual file name
nexus_data = Phylo.NexusIO.parse('D1_3_11-01-2021--07-04-2025_concat_318_tree_nexus.nex')

# Assuming there's one tree in the Nexus file
tree = next(nexus_data)

# Create a mapping of old labels to new labels
label_mapping = {
    # "OldLabel1": "NewLabelA",
    # "OldLabel2": "NewLabelB",
    # Add more mappings as needed
}

# Create new labels
os.chdir("C:/Users/maksiaevai.NCBI_NT/Documents/Avian_Flu_Files/")
county_info = pd.read_csv("Positive Case Info = Summary Counties.csv")

print(county_info["County"])
print(county_info["Event ID"])

for id in county_info["Event ID"].values:
    if id == id: # If not nan
        label_mapping[id] = county_info[county_info["Event ID"] == id]["County"].iloc[0] + "|" + county_info[county_info["Event ID"] == id]["Site Owner"].iloc[0].replace(" ", "_")


# Relabel the tips
for terminal in tree.get_terminals():
    # print(terminal.name)
    for key in label_mapping:
        if key in terminal.name:
            terminal.name = terminal.name.replace("'", "") + "|" + label_mapping[key] # Remove apostrophes
            print(terminal.name)

# Bootstrap values are usually associated with internal nodes
# Iterate through internal nodes and make sure their bootstrap values (node.comment or node.confidence) remain untouched during relabeling
for node in tree.get_nonterminals():
    # If bootstrap values are stored as node comments
    if node.comment and node.comment.startswith("bs"):
        # You might need to parse the comment to extract the value if it's not directly accessible
        # For simplicity, we assume they are already correctly stored as comments
        pass

    # If bootstrap values are stored as node confidence (e.g., node.confidence)
    if node.confidence:
        pass


# Save the relabeled tree to a new Nexus file
# Replace 'relabeled_tree.nex' with your desired output file name
Phylo.write([tree], 'D1_3_11-01-2021--07-04-2025_concat_318_tree_counties.nex', 'nexus')

