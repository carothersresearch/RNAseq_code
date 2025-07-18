import pandas as pd

def extract_info(attribute_str, key):
    """Extracts a value for a given key from GFF3 attribute column."""
    for entry in attribute_str.split(";"):
        if entry.startswith(key + "="):
            return entry.split("=")[1]
    return None

def extract_gene_lengths(gff_file, feature_type="CDS"):
    """
    Extracts gene lengths from a GFF3 file.
    :param gff_file: Path to the GFF3 file.
    :param feature_type: GFF feature type to use (default: CDS).
    :return: DataFrame with locus_tag and length.
    """
    # Read GFF, ignore comment lines
    df = pd.read_csv(gff_file, sep="\t", comment="#", header=None,
                     names=["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"])

    # Filter for feature type (e.g., CDS or gene)
    df = df[df["type"] == feature_type].copy()

    # Extract locus_tag
    df["locus_tag"] = df["attributes"].apply(lambda x: extract_info(x, "locus_tag"))

    # Calculate length
    df["length"] = df["end"] - df["start"] + 1

    # Group by locus_tag (some genes have multiple CDS parts)
    length_df = df.groupby("locus_tag")["length"].sum().reset_index()

    return length_df

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Extract gene lengths from GFF file.")
    parser.add_argument("gff_file", help="Input GFF3 file")
    parser.add_argument("-o", "--output", help="Output TSV file", default="gene_lengths.tsv")
    parser.add_argument("--feature", help="Feature type to use (CDS or gene)", default="CDS")
    args = parser.parse_args()

    result = extract_gene_lengths(args.gff_file, args.feature)
    result.to_csv(args.output, sep="\t", index=False)
    print(f"Gene lengths saved to {args.output}")
