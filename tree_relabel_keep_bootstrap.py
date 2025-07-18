from Bio import Phylo
from Bio import Nexus
from Bio import SeqIO

records = SeqIO.parse("THIS_IS_YOUR_INPUT_FILE.fasta", "fasta")
count = SeqIO.write(records, "THIS_IS_YOUR_OUTPUT_FILE.nexus", "nexus")

# Load the Nexus file
# Replace 'your_nexus_file.nex' with the actual file name
nexus_data = Phylo.NexusIO.parse('your_nexus_file.nex')

# Assuming there's one tree in the Nexus file
tree = next(nexus_data)

# Create a mapping of old labels to new labels
label_mapping = {
    "OldLabel1": "NewLabelA",
    "OldLabel2": "NewLabelB",
    # Add more mappings as needed
}

# Relabel the tips
for terminal in tree.get_terminals():
    if terminal.name in label_mapping:
        terminal.name = label_mapping[terminal.name]

# Bootstrap values are usually associated with internal nodes
# Iterate through internal nodes and make sure their bootstrap values (node.comment or node.confidence) remain untouched during relabeling
for node in tree.get_nonterminals():
    # If bootstrap values are stored as node comments
    if node.comment and node.comment.startswith("bootstrap"):
        # You might need to parse the comment to extract the value if it's not directly accessible
        # For simplicity, we assume they are already correctly stored as comments
        pass

    # If bootstrap values are stored as node confidence (e.g., node.confidence)
    if node.confidence:
        pass


# Save the relabeled tree to a new Nexus file
# Replace 'relabeled_tree.nex' with your desired output file name
Phylo.NexusIO.write(tree, 'relabeled_tree.nex', format='nexus')

