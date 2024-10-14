import os
from Bio import Entrez, SeqIO
from collections import defaultdict
import time

# Set your email
Entrez.email = "odw014@csiro.au"

#list of species to exclude
# All species have over 500 sequences on genbank and have complete genomes that are covered in refseq. They represent redundancy here and slow down the sorting/downloading process here
species_to_exclude = [
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
    "Bluetongue virus",
    "Monkeypox virus",
    "Orthoavulavirus javaense",
    "Avian coronavirus",
    "Protoparvovirus carnivoran1",
    "Human papillomavirus 16",
    "Bandavirus dabieense",
    "dengue virus type 3",
    "Torque teno midi virus",
    "Orthoflavivirus nilense",
    "Porcine epidemic diarrhea virus",
    "Betaarterivirus suid 2",
    "Chikungunya virus",
    "Zaire ebolavirus",
    "Foot-and-mouth disease virus",
    "Enterovirus C",
    "Enterovirus B",
    "human papillomavirus 31",
    "Adeno-associated virus",
    "Betacoronavirus 1",
    "human gammaherpesvirus 4",
    "Orthoflavivirus zikaense",
    "Enterovirus D",
    "Tomato yellow leaf curl virus",
    "Mammarenavirus lassaense",
    "Yellow fever virus",
    "Tomato leaf curl New Delhi virus",
    "Cytomegalovirus humanbeta5",
    "Human alphaherpesvirus 2",
    "Human mastadenovirus B",
    "Rhinovirus A",
    "Human alphaherpesvirus 1",
    "Orthobunyavirus oropoucheense",
    "Porcine circovirus 3",
    "Cucumber mosaic virus",
    "Polar freshwater circular DNA virus",
    "Orbivirus alphaequi",
    "Respirovirus pneumoniae",
    "Banana bunchy top virus",
    "Dengue virus",
    "Human mastadenovirus C",
    "Japanese encephalitis virus",
    "Mammalian orthoreovirus",
    "Human immunodeficiency virus 2",
    "Morbillivirus canis",
    "Isavirus salaris",
    "Metapneumovirus hominis",
    "African swine fever virus",
    "Piscine orthoreovirus",
    "Turnip mosaic virus",
    "Listeria phage LM 4-14-1",
    "Potato virus Y",
    "Simian foamy virus",
    "Porcine picobirnavirus",
    "Measles morbillivirus",
    "Pestivirus bovis",
    "Rabbit hemorrhagic disease virus",
    "Betapolyomavirus secuhominis",
    "Salmonella phage S.Kent 1-2-1",
    "Chicken anemia virus",
    "Orthoflavivirus encephalitidis",
    "Alphacoronavirus 1",
    "Alphapapillomavirus 7",
    "Infectious bursal disease virus",
    "Suid alphaherpesvirus 1",
    "Orthonairovirus haemorrhagiae",
    "Listeria phage LM 4-11-1",
    "Orthohantavirus puumalaense",
    # Additional species (previously added)
    "Mumps virus genotype G",
    "Orthotospovirus tomatomaculae",
    "Human papillomavirus",
    "Viral hemorrhagic septicemia virus",
    "Phlebovirus riftense",
    "Respiratory syncytial virus",
    "Orthohantavirus hantanense",
    "Middle East respiratory syndrome-related coronavirus",
    "Maize streak virus",
    "human papillomavirus 58",
    "Feline leukemia virus",
    "Rice black streaked dwarf virus",
    "human papillomavirus 35",
    "Influenza D virus",
    "Hepatitis delta virus",
    "Betapolyomavirus hominis",
    "Epizootic hemorrhagic disease virus",
    "Dependoparvovirus primate1",
    "human papillomavirus 52",
    "Sapporo virus",
    "Hepatitis E virus",
    "Orthohantavirus seoulense",
    "Usutu virus",
    "Bovine leukemia virus",
    "Equine infectious anemia virus",
    "Rhinovirus C",
    "Grapevine rupestris stem pitting-associated virus",
    "Plum pox virus",
    "Pigeon circovirus",
    "Respirovirus laryngotracheitidis",
    "Grapevine fanleaf virus",
    "Eastern equine encephalitis virus",
    "Primate T-lymphotropic virus 1",
    "Human mastadenovirus D",
    "Fort Morgan virus",
    "Infectious pancreatic necrosis virus",
    "Avian leukosis virus",
    "Alphapapillomavirus 10",
    "Iflavirus aladeformis",
    "Feline immunodeficiency virus",
    "Feline calicivirus",
    "Pepino mosaic virus",
    "Mumps orthorubulavirus",
    "Gallid alphaherpesvirus 1",
    "Enterovirus A",
    "dengue virus type 2",
    "TTV-like mini virus",
    "Enterobacteria phage SP",
    "uncultured Caudovirales phage",
    "Torque teno virus",
    "dengue virus type 4",
    "Enterovirus A",
    "dengue virus type 2",
    "TTV-like mini virus",
    "Enterobacteria phage SP",
    "uncultured Caudovirales phage",
    "Torque teno virus",
    "dengue virus type 4",
    "Avian orthoreovirus",
    "Rotavirus C",
    "Porcine circovirus 2",
    "human papillomavirus 35",
    "Rice black streaked dwarf virus",
    "unidentified influenza virus",
    "human papillomavirus 58",
    "uncultured human fecal virus",
    "uncultured Mediterranean phage uvMED",
    "Parechovirus A",
    "Cotton leaf curl Multan betasatellite",
    "uncultured marine virus",
    "Mogiana tick virus",
    "Emaravirus rosae",
    "Paslahepevirus balayani",
    "Staphylococcus phage SA 1298 Kay",
    "Begomovirus manihotis",
    "Rotavirus B",
    "Pegivirus hominis",
    "Human gammaherpesvirus 8",
    "Ovine progressive pneumonia virus",
    "Pestivirus suis",
    "Human alphaherpesvirus 3",
    "Beet necrotic yellow vein virus",
    "Infectious hematopoietic necrosis virus",
    "Ross River virus",
    "Beak and feather disease virus",
    "Southern rice black-streaked dwarf virus",
    "Potato virus X",
    "Fowl aviadenovirus C",
    "Equid alphaherpesvirus 1",
    "Human coronavirus 229E",
    "Chilli leaf curl virus",
    "Amdoparvovirus carnivoran1",
    "Protoparvovirus ungulate1",
    "Human coronavirus NL63",
    "Human mastadenovirus F",
    "Powassan virus",
    "Vaccinia virus",
    "Grapevine Pinot gris virus",
    "Bocaparvovirus primate1",
    "Erythroparvovirus primate1",
    "Hepatovirus A",
    "Chilli leaf curl Pakistan alphasatellite",
    "Dependoparvovirus anseriform1",
    "Tomato severe rugose virus",
    "Ageratum yellow vein Singapore alphasatellite",
    "Uukuvirus dabieshanense",
    "St. Louis encephalitis virus",
    "Human mastadenovirus E",
    "Duck circovirus",
    "Alfamovirus AMV",
    "Morbillivirus caprinae",
    "Feline parvovirus",
    "Human bocavirus",
    "Influenza C virus",
    "Coronavirus HKU15",
    "Mungbean yellow mosaic India virus",
    "Respirovirus bovis",
    "Rubella virus"]


# Construct the species exclusion query
species_exclusion_query = "NOT (" + " OR ".join([f'"{species}"[porgn]' for species in species_to_exclude]) + ")"

# Updated search query with length max 500,000 and species exclusions
search_query = f'(viruses[filter] AND is_nuccore[filter] AND ("1400"[SLEN] : "500000"[SLEN]) AND biomol_genomic[PROP]) {species_exclusion_query}'

# Perform the search
handle = Entrez.esearch(db="nuccore", term=search_query, retmax=500000)
record = Entrez.read(handle)
handle.close()

seq_ids = record['IdList']
total_sequences = len(seq_ids)
print(f"Total sequences retrieved: {total_sequences}")

# Function to batch process sequence IDs
def batch(iterable, n=500):
    """Yield successive n-sized batches from iterable."""
    l = len(iterable)
    for ndx in range(0, l, n):
        yield iterable[ndx:min(ndx + n, l)]

# Dictionary to hold TaxIDs and their sequence IDs
taxid_seq_ids = defaultdict(list)

# Batch size for fetching
batch_size = 500

# Counter for tracking processed sequences
processed_sequences = 0

# Fetch TaxIDs using ESummary
for seq_id_batch in batch(seq_ids, batch_size):
    ids = ','.join(seq_id_batch)
    handle = Entrez.esummary(db="nuccore", id=ids)
    records = Entrez.read(handle)
    handle.close()
    
    for record in records:
        seq_id = record['Id']
        taxid = record['TaxId']
        taxid_seq_ids[taxid].append(seq_id)
        processed_sequences += 1

    # Print a message every 5000 sequences
    if processed_sequences % 5000 == 0 or processed_sequences == total_sequences:
        print(f"Processed {processed_sequences} sequences out of {total_sequences}")

    # Respect NCBI rate limits
    time.sleep(0.05)

print(f"Total unique TaxIDs (species): {len(taxid_seq_ids)}")

# Select up to 50 sequences per species
selected_seq_ids = []

for taxid, ids in taxid_seq_ids.items():
    selected_ids = ids[:40]  # Get up to 40 IDs per species
    selected_seq_ids.extend(selected_ids)

print(f"Total selected sequences: {len(selected_seq_ids)}")

# Function to write sequences to a file
def download_sequences(seq_ids, filename, batch_size=500):
    processed = 0
    total = len(seq_ids)
    with open(filename, 'w') as out_handle:
        for seq_id_batch in batch(seq_ids, batch_size):
            ids = ','.join(seq_id_batch)
            handle = Entrez.efetch(db="nuccore", id=ids, rettype="fasta", retmode="text")
            data = handle.read()
            handle.close()
            out_handle.write(data)
            processed += len(seq_id_batch)

            # Print a message every 1000 sequences
            if processed % 1000 == 0 or processed == total:
                print(f"Downloaded {processed} sequences out of {total}")

            # Respect NCBI rate limits
            time.sleep(0.05)

print("Downloading selected sequences...")
output_file = "viral_sequences_selected.fasta"
download_sequences(selected_seq_ids, output_file)
print(f"Sequences saved to {output_file}")

# Function to extract accession numbers from sequence IDs
def get_accession_number(seq_id):
    """
    Extracts the accession number from a sequence ID.
    """
    # Split by dot to remove version number, if present
    return seq_id.split('.')[0]

# Function to build accession to TaxID mapping
def build_accession_to_taxid_mapping(input_fasta, mapping_files):
    # Extract accession numbers from the input_fasta
    accession_set = set()
    for record in SeqIO.parse(input_fasta, 'fasta'):
        seq_id = record.id
        accession = get_accession_number(seq_id)
        accession_set.add(accession)

    accession_to_taxid = {}

    for mapping_file in mapping_files:
        print(f"Processing mapping file: {mapping_file}")
        with open(mapping_file, 'r') as f:
            # Skip the header line if present
            header = f.readline()
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 3:
                    continue
                acc = parts[0]
                taxid = parts[2]
                if acc in accession_set:
                    accession_to_taxid[acc] = taxid
                    # Remove from accession_set to speed up
                    accession_set.remove(acc)
                    if not accession_set:
                        break
        if not accession_set:
            break

    if accession_set:
        print(f"Warning: TaxID not found for {len(accession_set)} accession(s)")

    return accession_to_taxid

working_directory = os.getcwd()

# Construct the paths to the mapping files
mapping_files = [
    os.path.join(working_directory, 'taxonomy', 'nucl_gb.accession2taxid'),
    os.path.join(working_directory, 'taxonomy', 'nucl_wgs.accession2taxid')
]
# Build the accession to TaxID mapping
accession_to_taxid = build_accession_to_taxid_mapping(output_file, mapping_files)

# Function to format FASTA headers for Kraken2
def format_headers_for_kraken(input_fasta, output_fasta, accession_to_taxid):
    processed = 0
    total = sum(1 for _ in SeqIO.parse(input_fasta, 'fasta'))
    with open(output_fasta, 'w') as out_handle:
        for record in SeqIO.parse(input_fasta, 'fasta'):
            seq_id = record.id
            accession = get_accession_number(seq_id)
            taxid = accession_to_taxid.get(accession, '0')  # Default to '0' if taxid not found
            # Format header
            record.id = f"{seq_id}|kraken:taxid|{taxid}"
            record.description = ''
            SeqIO.write(record, out_handle, 'fasta')
            processed += 1

            # Print a message every 1000 sequences
            if processed % 1000 == 0 or processed == total:
                print(f"Formatted {processed} sequences out of {total}")

    print(f"Formatted sequences saved to {output_fasta}")

# Format the headers
formatted_fasta = "viral_sequences_kraken.fasta"
format_headers_for_kraken(output_file, formatted_fasta, accession_to_taxid)