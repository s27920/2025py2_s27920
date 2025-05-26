from Bio import Entrez, SeqIO
from io import StringIO
import csv
import matplotlib.pyplot as plt

class NCBIRetriever:
    def __init__(self, email, api_key):
        Entrez.email = email
        Entrez.api_key = api_key
        Entrez.tool = 'BioScriptEx10'
    def search_taxid(self, taxid, min_len=None, max_len=None):
        print(f"Searching for records with taxID: {taxid}")
        try:
            handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
            records = Entrez.read(handle)
            organism_name = records[0]["ScientificName"]
            print(f"Organism: {organism_name} (TaxID: {taxid})")

            search_term = f"txid{taxid}[Organism] AND {min_len if min_len is not None else 0}:{max_len if max_len is not None else 10000000000000}[SLEN]"

            handle = Entrez.esearch(db="nucleotide", term=search_term, usehistory="y")
            search_results = Entrez.read(handle)
            count = int(search_results["Count"])

            if count == 0:
                print("No records found matching the criteria")
                return None

            print(f"Found {count} records")

            self.webenv = search_results["WebEnv"]
            self.query_key = search_results["QueryKey"]

            return count

        except Exception as e:
            print(f"Error searching TaxID {taxid}: {e}")
            return None

    def fetch_records(self, start=0, max_records=10):
        try:
            raw_text = Entrez.efetch(
                db="nucleotide",
                rettype="gb",
                retmode="text",
                retstart=start,
                retmax=max_records,
                webenv=self.webenv,
                query_key=self.query_key
            ).read()
            records = list(SeqIO.parse(StringIO(raw_text), "genbank"))
            return raw_text, records
        except Exception as e:
            print(f"Error fetching records: {e}")
            return "", []

def generate_csv_report(records, filename):
    with open(filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        for record in records:
            writer.writerow([record.annotations.get('accessions', [''])[0], len(record.seq), record.description])

def plot_sequence_lengths(records, filename):
    if not records:
        print("No records to plot")
        return

    sorted_records = sorted(records, key=lambda x: len(x.seq))
    accessions = [rec.annotations.get('accessions', [''])[0] for rec in sorted_records]

    plt.figure(figsize=(12, 6))
    plt.plot(range(len(accessions)), [len(rec.seq) for rec in sorted_records], marker='o', linestyle='-', color='b')
    plt.xticks(range(len(accessions)), accessions, rotation=45, ha='right')
    plt.ylabel('Sequence Length')
    plt.title('Sequence Length Distribution')
    plt.savefig(filename)
    plt.close()

def get_filters():
    min_len = None
    max_len = None
    try:
        min_input = input("Minimum sequence length [optional]: ").strip()
        max_input = input("Maximum sequence length [optional]: ").strip()
        min_len = int(min_input) if min_input else None
        max_len = int(max_input) if max_input else None
        if min_len and max_len and min_len > max_len:
            min_len, max_len = max_len, min_len
    except ValueError:
        print("Invalid filter")
    return min_len, max_len

def main():
    email = input("Enter your NCBI email: ")
    api_key = input("Enter your NCBI API key: ")

    retriever = NCBIRetriever(email, api_key)

    taxid = input("Enter taxonomic ID: ").strip()

    min_len, max_len = get_filters()

    count = retriever.search_taxid(taxid, min_len, max_len)

    raw_text, records = retriever.fetch_records(max_records=min(count, 500))
    if not records:
        print("No records fetched")
        return

    output_gb = f"taxid_{taxid}_records.gb"
    with open(output_gb, 'w') as f:
        f.write(raw_text)
    print(f"Saved {len(records)} records to {output_gb}")

    csv_file = f"taxid_{taxid}_report.csv"
    generate_csv_report(records, csv_file)
    print(f"Generated CSV report: {csv_file}")

    plot_file = f"taxid_{taxid}_lengths.png"
    plot_sequence_lengths(records, plot_file)
    print(f"Generated length plot: {plot_file}")

if __name__ == "__main__":
    main()