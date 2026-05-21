"""
find_closest_gene.py

For each indel in a BED file (chrom, start, end — 0-based half-open),
finds the closest downstream gene within 100,000 bp using a GFF file,
restricted to genes that appear in a featureCounts output file.

If the nearest downstream gene is not in the featureCounts list, it
counts how many genes are skipped before finding one that is.

Output TSV columns:
    indel_chrom, indel_start, indel_end, indel_length,
    gene_chrom, gene_start, gene_end, gene_name, genes_skipped
"""

import csv
from collections import defaultdict
from argparse import ArgumentParser

# ---------------------------------------------------------------------------
# 1. Parse the BED file  (0-based half-open: start inclusive, end exclusive)
# ---------------------------------------------------------------------------

def parse_bed(bed_path: str) -> list[tuple[str, int, int]]:
    """
    Read a BED file and return a list of (chrom, start, end) tuples.

    - BED coordinates are already 0-based half-open, so no conversion needed.
    - Lines beginning with '#', 'track', or 'browser' are skipped.
    - Only the first three columns are used; extra columns are ignored.
    """
    indels: list[tuple[str, int, int]] = []
    with open(bed_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith(("#", "track", "browser")):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            chrom = parts[0]
            start = int(parts[1])
            end   = int(parts[2])
            if start == end:
                length = int(parts[5]) - int(parts[4])
            else:
                length = int(parts[2]) - int(parts[1])
            indels.append((chrom, start, end, length))
    return indels


# ---------------------------------------------------------------------------
# 2. Parse the featureCounts output  →  set of gene names
# ---------------------------------------------------------------------------

def parse_aligned_genes(ag_path: str) -> set[str]:
    """
    Parse a featureCounts output file and return the set of gene IDs
    (Geneid column) present in it.

    featureCounts file structure
    ----------------------------
    Line 1  : '# Program:featureCounts ...'   ← comment, skip
    Line 2  : column headers (tab-separated), first column is 'Geneid'
    Line 3+ : data rows

    The gene ID is always in the first column (index 0).
    """
    gene_set: set[str] = set()
    with open(ag_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            parts = line.split("\t")
            # Skip the header row itself ("Geneid")
            if parts[0] == '"reference_gene_ID"':
                continue
            gene_id = parts[0].split(':')[1].strip('"')
            if gene_id:
                gene_set.add(gene_id)
    return gene_set


# ---------------------------------------------------------------------------
# 3. Parse the GFF file  →  genes sorted by position per chromosome
# ---------------------------------------------------------------------------

def parse_gff(gff_path: str) -> dict[str, list[dict]]:
    """
    Parse a GFF3 (or GTF) file and return a dict mapping
    chromosome -> list of gene records sorted by start position.

    Each record: {"chrom": str, "start": int, "end": int, "name": str}

    Only rows where column 3 == 'gene' are kept.
    GFF coordinates (1-based inclusive) are converted to 0-based half-open:
        start = col4 - 1,  end = col5  (unchanged)
    """
    genes_by_chrom: dict[str, list[dict]] = defaultdict(list)
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, _, feature, start_s, end_s, _, orientation, _, attrs = parts[:9]
            if feature != "gene":
                continue
            start = int(start_s) - 1   # GFF 1-based → 0-based
            end   = int(end_s)          # GFF inclusive end → half-open end
            name = _extract_gene_name(attrs)
            genes_by_chrom[chrom].append(
                {"chrom": chrom, "start": start, "end": end, "name": name, "orientation":orientation}
            )
    for chrom in genes_by_chrom:
        genes_by_chrom[chrom].sort(key=lambda g: g["start"])
    return dict(genes_by_chrom)


def _extract_gene_name(attrs: str) -> str:
    """
    Pull the gene name/ID out of a GFF3 or GTF attribute string.

    Priority order:
      GFF3 : Name=  >  ID=  (strips common 'gene:' prefix)
      GTF  : gene_id "..."
    """
    # GFF3
    for token in attrs.split(";"):
        token = token.strip()
        if token.startswith("Name="):
            return token[5:]
        if token.startswith("ID="):
            val = token[3:]
            return val.split(":")[-1]   # drop "gene:" prefix if present
    # GTF fallback
    for token in attrs.split(";"):
        token = token.strip()
        if token.startswith("gene_id"):
            return token.split('"')[1] if '"' in token else token.split()[-1]
    return attrs  # last resort


# ---------------------------------------------------------------------------
# 4. Core lookup: closest downstream gene that is in the featureCounts set
# ---------------------------------------------------------------------------

def find_closest_downstream_gene(
    indel_chrom: str,
    indel_start: int,
    indel_end: int,
    genes_by_chrom: dict[str, list[dict]],
    gene_set: set[str],
    max_distance: int = 100_000,
) -> dict | None:
    """
    Find the closest gene that:
      - is on the same chromosome as the indel
      - starts strictly after indel_end  (downstream only)
      - has its start within max_distance bp of indel_end
      - is present in gene_set
 
    Genes between the indel and the first matching gene are counted as
    'genes_skipped' (0 means the very first downstream gene matched).
 
    gene_location values:
        'within'     - the indel is entirely within the gene body (always returned
                       immediately, no other gene can win)
        'upstream'   - TSS is before the indel (strand-aware)
        'downstream' - TSS is after the indel (strand-aware)
 
    distance_to_TSS is always the distance from the indel start to the TSS,
    strand-aware, even for 'within' genes.
 
    Returns a result dict or None if nothing qualifies.
    """
    genes = genes_by_chrom.get(indel_chrom, [])
    if not genes:
        return None
 
    # Binary search: first gene whose start > indel_end
    lo, hi = 0, len(genes)
    while lo < hi:
        mid = (lo + hi) // 2
        if genes[mid]["start"] <= indel_end:
            lo = mid + 1
        else:
            hi = mid
    first_idx = lo
 
    genes_skipped = 0
    closest_downstream_gene = None
    for gene in genes[first_idx:]:
        if gene["start"] - indel_end > max_distance:
            break
        if gene["name"] in gene_set:
            closest_downstream_gene = {'chrom': gene["chrom"], 'start': gene["start"], 'end': gene["end"], 'name': gene["name"], 'orientation': gene["orientation"], 'genes_skipped': genes_skipped}
            break
        genes_skipped += 1
 
    closest_upstream_gene = None
    genes_skipped = 0
    for gene in genes[(first_idx - 1)::-1]:
        if gene['start'] < indel_start:
            if indel_start - gene["end"] > max_distance:
                break
            if gene["name"] in gene_set:
                closest_upstream_gene = {'chrom': gene["chrom"], 'start': gene["start"], 'end': gene["end"], 'name': gene["name"], 'orientation': gene["orientation"], 'genes_skipped': genes_skipped}
                break
            genes_skipped += 1
 
    if closest_upstream_gene and closest_downstream_gene:
        # check if indel is within the upstream gene body — return immediately if so
        if indel_start >= closest_upstream_gene['start'] and indel_end <= closest_upstream_gene['end']:
            if closest_upstream_gene['orientation'] == '+':
                upstream_start_distance = indel_start - closest_upstream_gene['start']
            else:
                upstream_start_distance = closest_upstream_gene['end'] - indel_end
            return {
                "gene_chrom":                closest_upstream_gene["chrom"],
                "gene_start":                closest_upstream_gene["start"],
                "gene_end":                  closest_upstream_gene["end"],
                "gene_name":                 closest_upstream_gene["name"],
                "gene_orientation":          closest_upstream_gene["orientation"],
                "genes_skipped":             closest_upstream_gene["genes_skipped"],
                "distance_to_TSS_start_site": upstream_start_distance,
                "gene_location":             'within',
            }
        elif closest_upstream_gene['orientation'] == '+':
            upstream_gene_location = 'upstream'
            upstream_start_distance = indel_start - closest_upstream_gene['start']
        else:
            upstream_gene_location = 'upstream'
            upstream_start_distance = indel_start - closest_upstream_gene['end']
 
        # check if indel is within the downstream gene body — return immediately if so
        if indel_start >= closest_downstream_gene['start'] and indel_end <= closest_downstream_gene['end']:
            if closest_downstream_gene['orientation'] == '+':
                downstream_start_distance = indel_start - closest_downstream_gene['start']
            else:
                downstream_start_distance = closest_downstream_gene['end'] - indel_end
            return {
                "gene_chrom":                closest_downstream_gene["chrom"],
                "gene_start":                closest_downstream_gene["start"],
                "gene_end":                  closest_downstream_gene["end"],
                "gene_name":                 closest_downstream_gene["name"],
                "gene_orientation":          closest_downstream_gene["orientation"],
                "genes_skipped":             closest_downstream_gene["genes_skipped"],
                "distance_to_TSS_start_site": downstream_start_distance,
                "gene_location":             'within',
            }
        elif closest_downstream_gene['orientation'] == '+':
            downstream_gene_location = 'downstream'
            downstream_start_distance = closest_downstream_gene['start'] - indel_end
        else:
            downstream_gene_location = 'downstream'
            downstream_start_distance = closest_downstream_gene['end'] - indel_end
 
        if upstream_start_distance == min(upstream_start_distance, downstream_start_distance):
            return {
                "gene_chrom":                closest_upstream_gene["chrom"],
                "gene_start":                closest_upstream_gene["start"],
                "gene_end":                  closest_upstream_gene["end"],
                "gene_name":                 closest_upstream_gene["name"],
                "gene_orientation":          closest_upstream_gene["orientation"],
                "genes_skipped":             closest_upstream_gene["genes_skipped"],
                "distance_to_TSS_start_site": upstream_start_distance,
                "gene_location":             upstream_gene_location,
            }
        else:
            return {
                "gene_chrom":                closest_downstream_gene["chrom"],
                "gene_start":                closest_downstream_gene["start"],
                "gene_end":                  closest_downstream_gene["end"],
                "gene_name":                 closest_downstream_gene["name"],
                "gene_orientation":          closest_downstream_gene["orientation"],
                "genes_skipped":             closest_downstream_gene["genes_skipped"],
                "distance_to_TSS_start_site": downstream_start_distance,
                "gene_location":             downstream_gene_location,
            }
 
    elif closest_upstream_gene:
        # check if indel is within the upstream gene body — return immediately if so
        if indel_start >= closest_upstream_gene['start'] and indel_end <= closest_upstream_gene['end']:
            if closest_upstream_gene['orientation'] == '+':
                upstream_start_distance = indel_start - closest_upstream_gene['start']
            else:
                upstream_start_distance = closest_upstream_gene['end'] - indel_end
            return {
                "gene_chrom":                closest_upstream_gene["chrom"],
                "gene_start":                closest_upstream_gene["start"],
                "gene_end":                  closest_upstream_gene["end"],
                "gene_name":                 closest_upstream_gene["name"],
                "gene_orientation":          closest_upstream_gene["orientation"],
                "genes_skipped":             closest_upstream_gene["genes_skipped"],
                "distance_to_TSS_start_site": upstream_start_distance,
                "gene_location":             'within',
            }
        elif closest_upstream_gene['orientation'] == '+':
            upstream_gene_location = 'upstream'
            upstream_start_distance = indel_start - closest_upstream_gene['start']
        else:
            upstream_gene_location = 'upstream'
            upstream_start_distance = indel_start - closest_upstream_gene['end']
 
        return {
            "gene_chrom":                closest_upstream_gene["chrom"],
            "gene_start":                closest_upstream_gene["start"],
            "gene_end":                  closest_upstream_gene["end"],
            "gene_name":                 closest_upstream_gene["name"],
            "gene_orientation":          closest_upstream_gene["orientation"],
            "genes_skipped":             closest_upstream_gene["genes_skipped"],
            "distance_to_TSS_start_site": upstream_start_distance,
            "gene_location":             upstream_gene_location,
        }
 
    elif closest_downstream_gene:
        # check if indel is within the downstream gene body — return immediately if so
        if indel_start >= closest_downstream_gene['start'] and indel_end <= closest_downstream_gene['end']:
            if closest_downstream_gene['orientation'] == '+':
                downstream_start_distance = indel_start - closest_downstream_gene['start']
            else:
                downstream_start_distance = closest_downstream_gene['end'] - indel_end
            return {
                "gene_chrom":                closest_downstream_gene["chrom"],
                "gene_start":                closest_downstream_gene["start"],
                "gene_end":                  closest_downstream_gene["end"],
                "gene_name":                 closest_downstream_gene["name"],
                "gene_orientation":          closest_downstream_gene["orientation"],
                "genes_skipped":             closest_downstream_gene["genes_skipped"],
                "distance_to_TSS_start_site": downstream_start_distance,
                "gene_location":             'within',
            }
        elif closest_downstream_gene['orientation'] == '+':
            downstream_gene_location = 'downstream'
            downstream_start_distance = closest_downstream_gene['start'] - indel_end
        else:
            downstream_gene_location = 'downstream'
            downstream_start_distance = closest_downstream_gene['end'] - indel_end
 
        return {
            "gene_chrom":                closest_downstream_gene["chrom"],
            "gene_start":                closest_downstream_gene["start"],
            "gene_end":                  closest_downstream_gene["end"],
            "gene_name":                 closest_downstream_gene["name"],
            "gene_orientation":          closest_downstream_gene["orientation"],
            "genes_skipped":             closest_downstream_gene["genes_skipped"],
            "distance_to_TSS_start_site": downstream_start_distance,
            "gene_location":             downstream_gene_location,
        }
 
    return None

# ---------------------------------------------------------------------------
# 5. Process all indels
# ---------------------------------------------------------------------------

def process_indels(
    indels: list[tuple[str, int, int]],
    gff_path: str,
    aligned_genes_path: str,
    max_distance: int = 100_000,
) -> list[dict]:
    """
    Parameters
    ----------
    indels            : list of (chrom, start, end) from parse_bed()
    gff_path          : path to GFF/GFF3 annotation file
    featurecounts_path: path to featureCounts output (gene whitelist)
    max_distance      : upstream search limit in bp (default 100,000)

    Returns
    -------
    List of result dicts (one per indel).  Indels with no match have
    None in all gene_* fields and genes_skipped.
    """
    genes_by_chrom = parse_gff(gff_path)
    gene_set       = parse_aligned_genes(aligned_genes_path)
    print(f"Loaded {sum(len(v) for v in genes_by_chrom.values())} genes "
          f"across {len(genes_by_chrom)} chromosomes from GFF.")
    print(f"Loaded {len(gene_set)} gene IDs from aligned gene file.")
    results = []
    for chrom, start, end, indel_length in indels:
        hit = find_closest_downstream_gene(
            chrom, start, end, genes_by_chrom, gene_set, max_distance
        )
        results.append({
            "indel_chrom":   chrom,
            "indel_start":   start,
            "indel_end":     end,
            "indel_length":  indel_length,
            "gene_chrom":    hit["gene_chrom"]    if hit else None,
            "gene_start":    hit["gene_start"]    if hit else None,
            "gene_end":      hit["gene_end"]      if hit else None,
            "gene_name":     hit["gene_name"]     if hit else None,
            "gene_orientation":     hit["gene_orientation"]     if hit else None,
            "genes_skipped": hit["genes_skipped"] if hit else None,
            "distance_to_TSS": hit["distance_to_TSS_start_site"] if hit else None,
            "gene_location": hit["gene_location"] if hit else None,
        })
    return results



# ---------------------------------------------------------------------------
# 6. Write TSV output
# ---------------------------------------------------------------------------

FIELDNAMES = [
    "indel_chrom", "indel_start", "indel_end", "indel_length",
    "gene_chrom",  "gene_start",  "gene_end",  "gene_name", "gene_orientation",
    "genes_skipped", "distance_to_TSS", "gene_location",
]

def write_results(results: list[dict], output_path: str) -> None:
    with open(output_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDNAMES, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(results)
    matched = sum(1 for r in results if r["gene_name"] is not None)
    print(f"Wrote {len(results)} rows ({matched} matched, "
          f"{len(results) - matched} unmatched) to {output_path}")


# ---------------------------------------------------------------------------
# 7. Entry point
# ---------------------------------------------------------------------------
parser = ArgumentParser()
#inputs
parser.add_argument("-indel_double_bed", "--indel_double_bed", help="File containing double_bed_coordinates of all indels", required=True)
parser.add_argument("-GFF", "--GFF", help="The reference gff file", required=True)
parser.add_argument("-aligned_genes", "--aligned_genes", help="File containing the genes that i could aligned between the refernece and alterante with RNA reads associated", required=True)

#outputs
parser.add_argument("-Output_file", "--Output_file", action='store', help="For each indel what is its closest gene", required=False)

args = parser.parse_args()



if __name__ == "__main__":
    # ── Edit these three paths ──────────────────────────────────────────────
    BED_FILE           = args.indel_double_bed
    GFF_FILE           = args.GFF
    ALIGNED_GENES_FILE = args.aligned_genes
    OUTPUT_FILE        = args.Output_file
    # ────────────────────────────────────────────────────────────────────────

    indels  = parse_bed(BED_FILE)
    print(f"Loaded {len(indels)} indels from {BED_FILE}.")

    results = process_indels(indels, GFF_FILE, ALIGNED_GENES_FILE)
    write_results(results, OUTPUT_FILE)
