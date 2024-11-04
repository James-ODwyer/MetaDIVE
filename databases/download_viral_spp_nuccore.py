import os
from Bio import Entrez, SeqIO
from collections import defaultdict
import time
import datetime

#set email
Entrez.email = "odw014@csiro.au"

#list of species to exclude
# All species have over 500 sequences on genbank and have complete genomes that are covered in refseq. They represent redundancy here and slow down the sorting/downloading process here
species_to_exclude = [
    "Protoparvovirus carnivoran1", "Human papillomavirus 16", "avian paramyxovirus 1",
    "Enterovirus A", "Avian coronavirus", "dengue virus type 1", "Bluetongue virus",
    "Lyssavirus rabies", "Bacteriophage sp.", "Norwalk virus", "Simian-Human immunodeficiency virus",
    "Hepatitis B virus", "Rotavirus A", "Simian immunodeficiency virus", "Influenza B virus",
    "Caudoviricetes sp.", "Hepacivirus hominis", "Influenza A virus", "Human immunodeficiency virus 1",
    "Severe acute respiratory syndrome-related coronavirus", "Circoviridae sp.",
    "uncultured Caudovirales phage", "Riboviria sp.", "dengue virus type 3",
    "Foot-and-mouth disease virus", "Zaire ebolavirus", "Chikungunya virus", "Enterovirus C",
    "Betaarterivirus suid 2", "Torque teno virus", "Porcine epidemic diarrhea virus", "Porcine circovirus 2",
    "Orthoflavivirus nilense", "Torque teno midi virus", "Bandavirus dabieense", "TTV-like mini virus",
    "dengue virus type 2", "Microviridae sp.", "Monkeypox virus", "Orthopneumovirus hominis",
    "Tomato yellow leaf curl virus", "Polar freshwater circular DNA virus", "Mammarenavirus lassaense",
    "Human mastadenovirus C", "Cytomegalovirus humanbeta5", "Orbivirus alphaequi", "Morbillivirus canis",
    "Leviviridae sp.", "Human mastadenovirus B", "Orthoflavivirus zikaense", "Enterovirus D",
    "Betacoronavirus 1", "Respirovirus pneumoniae", "Japanese encephalitis virus", "human gammaherpesvirus 4",
    "Adeno-associated virus", "Human alphaherpesvirus 1", "human papillomavirus 31",
    "Human immunodeficiency virus", "Enterovirus B", "Porcine picobirnavirus", "Porcine circovirus 3",
    "Orthonairovirus haemorrhagiae", "Orthobunyavirus oropoucheense", "Cucumber mosaic virus",
    "Mumps virus genotype G", "Alphainfluenzavirus influenzae", "Microvirus sp.", "Listeria phage LM 4-14-1",
    "uncultured virus", "Rhinovirus A", "dengue virus type 4", "Tomato leaf curl New Delhi virus",
    "Cressdnaviricota sp.", "Rabbit hemorrhagic disease virus", "Human alphaherpesvirus 2",
    "Measles morbillivirus", "Metapneumovirus hominis", "Yellow fever virus", "Avian orthoreovirus",
    "Alphapapillomavirus 7", "Maize streak virus", "Piscine orthoreovirus", "Anelloviridae sp.",
    "Orthohantavirus hantanense", "Middle East respiratory syndrome-related coronavirus",
    "Listeria phage LM 4-11-1", "Betapolyomavirus hominis", "uncultured Mediterranean phage uvMED",
    "Potato virus Y", "Human papillomavirus", "Dengue virus", "Picornavirales sp.",
    "Respiratory syncytial virus", "Rotavirus C", "Mammalian orthoreovirus", "Turnip mosaic virus",
    "Viral hemorrhagic septicemia virus", "Phlebovirus riftense", "Human immunodeficiency virus 2",
    "Rhinovirus C", "Grapevine rupestris stem pitting-associated virus", "Pigeon circovirus",
    "Chicken anemia virus", "human papillomavirus 52", "Sapporo virus", "Suid alphaherpesvirus 1",
    "Feline leukemia virus", "Salmonella phage S.Kent 1-2-1", "Dependoparvovirus primate1",
    "African swine fever virus", "Influenza D virus", "Epizootic hemorrhagic disease virus", "Mitovirus sp.",
    "uncultured human fecal virus", "Hepatitis delta virus", "human papillomavirus 58",
    "Orthotospovirus tomatomaculae", "human papillomavirus 35", "Rice black streaked dwarf virus",
    "Escherichia phage ECO95 1-6-1", "Alphapapillomavirus 10", "Feline calicivirus",
    "Orthoflavivirus encephalitidis", "Picobirnavirus sp.", "Plum pox virus", "Bovine leukemia virus",
    "Orthohantavirus puumalaense", "Avian leukosis virus", "Eastern equine encephalitis virus",
    "Human mastadenovirus D", "Paslahepevirus balayani", "Usutu virus", "Grapevine fanleaf virus",
    "Orthohantavirus seoulense", "Hepatitis E virus", "Betapolyomavirus secuhominis", "Alphacoronavirus 1",
    "unidentified influenza virus", "Infectious bursal disease virus", "Dicistroviridae sp.",
    "Human gammaherpesvirus 8", "Picornaviridae sp.", "Reoviridae sp.", "Infectious hematopoietic necrosis virus",
    "Staphylococcus phage SA 1298 Kay", "Feline immunodeficiency virus", "Pepino mosaic virus",
    "Human alphaherpesvirus 3", "uncultured marine virus", "Inoviridae sp.", "Southern rice black-streaked dwarf virus",
    "Totiviridae sp.", "Beak and feather disease virus", "Infectious pancreatic necrosis virus",
    "Respirovirus laryngotracheitidis", "Mogiana tick virus", "Parechovirus A", "Mumps orthorubulavirus",
    "Genomoviridae sp.", "Amdoparvovirus carnivoran1", "Gallid alphaherpesvirus 1", "Chilli leaf curl virus",
    "Protoparvovirus ungulate1", "Powassan virus", "Human adenovirus sp.", "Ambivirus sp.", "Hepatovirus A",
    "Emaravirus rosae", "Human coronavirus 229E", "Ovine progressive pneumonia virus", "Erythroparvovirus primate1",
    "Circular genetic element sp.", "Pestivirus suis", "Pestivirus bovis", "Grapevine Pinot gris virus",
    "Narnaviridae sp.", "Bocaparvovirus primate1", "Begomovirus manihotis", "Primate T-lymphotropic virus 1",
    "Feline parvovirus", "Iflaviridae sp.", "Alfamovirus AMV", "Human bocavirus", "Sindbis virus"]


# second smaller exclusion list for searching genomes db
species_to_exclude2 = [
    # Existing species (already in the script)
    "Severe acute respiratory syndrome-related coronavirus",
    "Alphainfluenzavirus influenzae",
    "Human immunodeficiency virus 1",
    "Betainfluenzavirus influenzae",
    "Hepacivirus hominis",
    "Rotavirus A",
    "Simian immunodeficiency virus",
    "Hepatitis B virus",
    "Orthopneumovirus hominis",
    "Simian-Human immunodeficiency virus",
    "Human immunodeficiency virus",
    "dengue virus type 1",
    "Lyssavirus rabies",
    "Norwalk virus",
    "Bluetongue virus"]

working_directory = os.getcwd()
taxonomy_dir = os.path.join(working_directory, "taxonomy")
nodes_path = os.path.join(taxonomy_dir, "nodes.dmp")
virus_root_taxid = "10239"

# Step 1: Extract viral TaxIDs from nodes.dmp
def extract_viral_taxids(nodes_file, virus_root_taxid):
    viral_taxids = set()
    with open(nodes_file, 'r') as f:
        for line in f:
            parts = line.split("\t|\t")
            taxid = parts[0].strip()
            parent_taxid = parts[1].strip()
            # Check if TaxID is a descendant of the virus root
            if parent_taxid == virus_root_taxid or parent_taxid in viral_taxids:
                viral_taxids.add(taxid)
    return viral_taxids

# Extract viral TaxIDs
viral_taxids = extract_viral_taxids(nodes_path, virus_root_taxid)
print(f"Total viral TaxIDs found: {len(viral_taxids)}")


viral_taxids = list(viral_taxids)

# Set batch parameters

# Function to split iterable into chunks
def batch(iterable, n=200):
    l = len(iterable)
    for ndx in range(0, l, n):
        yield iterable[ndx:min(ndx + n, l)]


# Construct the species exclusion query using species names only
species_exclusion_query = "NOT (" + " OR ".join([f'"{species}"[Organism]' for species in species_to_exclude]) + ")"

# Function to download nucleotide viruses with dynamic batching strategy and species exclusion
def download_nucleotide_viruses(viral_taxids):
    all_seq_ids = []
    taxid_count = 0  # Counter for tracking TaxIDs processed
    max_retries = 3  # Maximum number of retries for errors

    # Convert viral_taxids to a set initially for efficient removal
    viral_taxids_set = set(viral_taxids)
    
    # Record the start time of the entire download process
    start_time = time.time()

    # Retmax thresholds and corresponding batch levels
    Retmax1 = 6000
    Retmax2 = 4000
    Retmax3 = 1500
    Retmax4 = 400
    Retmax5 = 60

    # Process TaxIDs starting with the largest batch size of 120
    while viral_taxids_set:
        taxid_batch = list(viral_taxids_set)[:120]  # Get the first 120 taxids
        taxid_query = " OR ".join([f"txid{taxid}[Organism]" for taxid in taxid_batch])
        complete_query = f"({taxid_query}) {species_exclusion_query}"
        
        try:
            handle = Entrez.esearch(db="nuccore", term=complete_query, retmax=Retmax1)
            record = Entrez.read(handle)
            handle.close()
            count = int(record["Count"])

            # If count exceeds Retmax1, reduce batch size and retmax progressively
            if count > Retmax1:
                print(f"Batch of 120 too large ({count} sequences). Splitting into smaller batches of 30 with retmax={Retmax2}.")
                
                for sub_batch in batch(taxid_batch, 30):
                    sub_query = " OR ".join([f"txid{taxid}[Organism]" for taxid in sub_batch])
                    complete_query = f"({sub_query}) {species_exclusion_query}"
                    
                    handle = Entrez.esearch(db="nuccore", term=complete_query, retmax=Retmax2)
                    sub_record = Entrez.read(handle)
                    handle.close()
                    sub_count = int(sub_record["Count"])

                    # If count still exceeds Retmax2, go to the next level
                    if sub_count > Retmax2:
                        print(f"Sub-batch of 30 too large ({sub_count} sequences). Splitting into smaller batches of 10 with retmax={Retmax3}.")
                        
                        for mini_batch in batch(sub_batch, 10):
                            mini_query = " OR ".join([f"txid{taxid}[Organism]" for taxid in mini_batch])
                            complete_query = f"({mini_query}) {species_exclusion_query}"
                            
                            handle = Entrez.esearch(db="nuccore", term=complete_query, retmax=Retmax3)
                            mini_record = Entrez.read(handle)
                            handle.close()
                            mini_count = int(mini_record["Count"])

                            # If count still exceeds Retmax3, go to the next level
                            if mini_count > Retmax3:
                                print(f"Mini-batch of 10 too large ({mini_count} sequences). Splitting into pairs with retmax={Retmax4}.")
                                
                                for pair_batch in batch(mini_batch, 2):
                                    pair_query = " OR ".join([f"txid{taxid}[Organism]" for taxid in pair_batch])
                                    complete_query = f"({pair_query}) {species_exclusion_query}"
                                    
                                    handle = Entrez.esearch(db="nuccore", term=complete_query, retmax=Retmax4)
                                    pair_record = Entrez.read(handle)
                                    handle.close()
                                    pair_count = int(pair_record["Count"])

                                    # If count still exceeds Retmax4, process each TaxID individually
                                    if pair_count > Retmax4:
                                        print(f"Pair too large ({pair_count} sequences). Processing each TaxID individually with retmax={Retmax5}.")
                                        
                                        for taxid in pair_batch:
                                            individual_query = f"txid{taxid}[Organism] {species_exclusion_query}"
                                            handle = Entrez.esearch(db="nuccore", term=individual_query, retmax=Retmax5)
                                            individual_record = Entrez.read(handle)
                                            handle.close()
                                            all_seq_ids.extend(individual_record["IdList"])
                                            print(f"Retrieved {len(individual_record['IdList'])} sequences for TaxID {taxid}")
                                    else:
                                        all_seq_ids.extend(pair_record["IdList"])
                                        print(f"Retrieved {len(pair_record['IdList'])} sequences for pair batch")
                            else:
                                all_seq_ids.extend(mini_record["IdList"])
                                print(f"Retrieved {len(mini_record['IdList'])} sequences for mini-batch")
                    else:
                        all_seq_ids.extend(sub_record["IdList"])
                        print(f"Retrieved {len(sub_record['IdList'])} sequences for sub-batch")
            else:
                # If initial batch is within Retmax1, add all IDs to results
                all_seq_ids.extend(record["IdList"])
                print(f"Retrieved {count} sequences for batch of 120 TaxIDs")

            # Update counter and report progress
            taxid_count += len(taxid_batch)
            if taxid_count % 1000 == 0:
                elapsed_time = time.time() - start_time
                formatted_time = str(datetime.timedelta(seconds=int(elapsed_time)))
                print(f"Processed {taxid_count} TaxIDs. Elapsed time: {formatted_time}")

        except Exception as e:
            print(f"Error retrieving batch: {e}")
            time.sleep(5)  # Delay before retrying the batch

        # Remove processed taxids from the set to avoid reprocessing
        viral_taxids_set.difference_update(taxid_batch)
        
        # Small delay to respect NCBI server rate limits
        time.sleep(0.05)

    return all_seq_ids

# Run the function to get all nucleotide sequence IDs using viral_taxids
all_seq_ids = download_nucleotide_viruses(viral_taxids)
print(f"Total unique sequence IDs retrieved: {len(all_seq_ids)}")


# Function to fetch summaries with retry logic
def fetch_summary_with_retries(seq_id_batch, retries=3):
    """Fetch Entrez summaries with retries."""
    for attempt in range(retries):
        try:
            ids = ','.join(seq_id_batch)
            handle = Entrez.esummary(db="nuccore", id=ids)
            records = Entrez.read(handle)
            handle.close()
            return records
        except Exception as e:
            print(f"Attempt {attempt + 1} failed: {e}")
            if attempt < retries - 1:
                print("Retrying...")
                time.sleep(2)  # Delay before retrying
            else:
                print("Max retries reached. Skipping this batch.")
                return None

# Process sequence IDs in batches to retrieve TaxIDs with retry logic
taxid_seq_ids = defaultdict(list)
for seq_id_batch in batch(all_seq_ids, 500):
    records = fetch_summary_with_retries(seq_id_batch)
    if records is None:
        continue

    for record in records:
        seq_id = record['Id']
        taxid = record['TaxId']
        taxid_seq_ids[taxid].append(seq_id)

print(f"Total unique TaxIDs (species): {len(taxid_seq_ids)}")

# Select up to 60 sequences per species
selected_seq_ids = []
for taxid, ids in taxid_seq_ids.items():
    selected_ids = ids[:60]
    selected_seq_ids.extend(selected_ids)

print(f"Total selected sequences: {len(selected_seq_ids)}")

# Function to download sequences
def download_sequences(seq_ids, filename, batch_size=500):
    with open(filename, 'w') as out_handle:
        for seq_id_batch in batch(seq_ids, batch_size):
            ids = ','.join(seq_id_batch)
            retries = 0
            while retries < 3:
                try:
                    handle = Entrez.efetch(db="nuccore", id=ids, rettype="fasta", retmode="text")
                    data = handle.read()
                    handle.close()
                    out_handle.write(data)
                    print(f"Downloaded batch of {len(seq_id_batch)} sequences")
                    break
                except Exception as e:
                    retries += 1
                    print(f"Error during download, retrying ({retries}/3): {e}")
                    time.sleep(1)
            if retries == 3:
                print("Max retries reached for this batch, skipping.")

# Download selected sequences to file
output_file = "viral_combined_sequences.fasta"
download_sequences(selected_seq_ids, output_file)
print(f"Sequences saved to {output_file}")

output_file = "viral_combined_sequences.fasta"

# Function to extract accession numbers from sequence IDs
def get_accession_number(seq_id):
    """Extracts the accession number from a sequence ID by splitting out the version if present."""
    return seq_id.split('.')[0]

# Define path to the taxonomy mapping files for accessions and taxids
working_directory = os.getcwd()
mapping_files = [
    os.path.join(working_directory, 'taxonomy', 'nucl_gb.accession2taxid'),
    os.path.join(working_directory, 'taxonomy', 'nucl_wgs.accession2taxid')
]

# Function to build accession-to-TaxID mapping from taxonomy files
def build_accession_to_taxid_mapping(input_fasta, mapping_files):
    accession_set = set()
    for record in SeqIO.parse(input_fasta, 'fasta'):
        seq_id = record.id
        accession = get_accession_number(seq_id)
        accession_set.add(accession)

    accession_to_taxid = {}
    
    for mapping_file in mapping_files:
        print(f"Processing mapping file: {mapping_file}")
        with open(mapping_file, 'r') as f:
            header = f.readline()  # Skip header
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 3:
                    continue
                acc = parts[0]
                taxid = parts[2]
                if acc in accession_set:
                    accession_to_taxid[acc] = taxid
                    accession_set.remove(acc)
                    if not accession_set:
                        break
        if not accession_set:
            break

    if accession_set:
        print(f"Warning: TaxID not found for {len(accession_set)} accession(s)")

    return accession_to_taxid

# Build the accession-to-TaxID mapping for the downloaded sequences
accession_to_taxid = build_accession_to_taxid_mapping(output_file, mapping_files)

# Function to format FASTA headers for Kraken2, saving output with TaxIDs
def format_headers_for_kraken(input_fasta, output_fasta, accession_to_taxid):
    processed = 0
    total = sum(1 for _ in SeqIO.parse(input_fasta, 'fasta'))
    
    with open(output_fasta, 'w') as out_handle:
        for record in SeqIO.parse(input_fasta, 'fasta'):
            seq_id = record.id
            accession = get_accession_number(seq_id)
            taxid = accession_to_taxid.get(accession, '0')  # Default to '0' if taxid not found
            record.id = f"{seq_id}|kraken:taxid|{taxid}"
            record.description = ''
            SeqIO.write(record, out_handle, 'fasta')
            processed += 1

            if processed % 1000 == 0 or processed == total:
                print(f"Formatted {processed} sequences out of {total}")

    print(f"Formatted sequences saved to {output_fasta}")

# Format the headers for Kraken2
formatted_fasta = "viral_combined_sequences_kraken.fasta"
format_headers_for_kraken(output_file, formatted_fasta, accession_to_taxid)
print(f"Formatted sequences saved to {formatted_fasta}")