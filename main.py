# ============================================
# SARS-COV-2 GENOMIC ANALYSIS PIPELINE
# A complete tool to compare COVID-19 variants
# ============================================

# ============================================
# STEP 1: IMPORT ALL THE TOOLS WE NEED
# ============================================
# These are like downloading apps for your phone - each one does something specific

from Bio import SeqIO  
# BioPython - helps us read genetic sequence files (like FASTA)

from Bio.Align import PairwiseAligner  
# This helps us compare two sequences side by side to find differences

import pandas as pd  
# Pandas = "Python Data" - helps us make nice tables (like Excel in Python)

import numpy as np  
# NumPy = "Number Python" - helps with math calculations

from tabulate import tabulate  
# Makes our tables look pretty when printed

import matplotlib.pyplot as plt  
# Matplotlib = basic chart maker - creates graphs and figures

import seaborn as sns  
# Seaborn = makes even prettier charts than matplotlib

import os  
# Helps us work with files and folders (like creating directories)

# Set up a nice style for our charts so they look professional
plt.style.use('seaborn-v0_8-darkgrid')  # This just makes charts look better
sns.set_palette("husl")  # Picks nice colors for our graphs automatically


# ============================================
# STEP 2: READ ALL OUR COVID SEQUENCE FILES
# ============================================
# This is like opening each file and copying the genetic code into Python

print("📂 Reading sequence files...")

# For each file, we:
# 1. Open the file
# 2. Read the sequence
# 3. Store it in a variable (like a labeled box) so we can use it later

for record in SeqIO.parse("COVID_Sequence(Reference).fasta", "fasta"):
    reference_sequence = str(record.seq)  # Store reference genome
    print("  ✅ Loaded Reference (Wuhan-Hu-1)")

for record2 in SeqIO.parse("COVID_Sequence(Alpha Variant).fasta", "fasta"):
    alpha_variant = str(record2.seq)  # Store Alpha variant
    print("  ✅ Loaded Alpha Variant")

for record3 in SeqIO.parse("COVID_delta Variant.fasta", "fasta"):
    delta_variant = str(record3.seq)  # Store Delta variant
    print("  ✅ Loaded Delta Variant")

for record4 in SeqIO.parse("Omicron_(BA.1)COVID_Variant.fasta", "fasta"):
    Omicron_BA_1 = str(record4.seq)  # Store Omicron BA.1
    print("  ✅ Loaded Omicron BA.1")

for record5 in SeqIO.parse("Omicron_(XBB)_COVID_Variant.fasta", "fasta"):
    Omicron_XBB = str(record5.seq)  # Store Omicron XBB
    print("  ✅ Loaded Omicron XBB")


# ============================================
# STEP 3: CALCULATE GC CONTENT FOR ANY SEQUENCE
# ============================================
# GC content = how many G and C letters in the DNA
# This matters because GC bonds are stronger than AT bonds
# Higher GC = more stable genome

def GC_content_calculator(sequence):
    """
    This function takes a DNA sequence and figures out what percentage is G and C.
    It's like counting how many red marbles in a bag of mixed colors.
    """
    
    # Convert to uppercase just to be safe (files sometimes have lowercase)
    seq = str(sequence).upper()
    
    # Count each letter
    g_count = seq.count("G")  # How many G's?
    c_count = seq.count("C")  # How many C's?
    a_count = seq.count("A")  # How many A's?
    t_count = seq.count("T")  # How many T's?
    
    # Valid bases = all the letters that are actually A,T,G,C (not N or other)
    valid_bases = a_count + t_count + g_count + c_count
    
    # Safety check - if there are no valid bases, return 0 (avoid division by zero)
    if valid_bases == 0:
        return 0
    
    # Calculate percentage: (G+C) divided by total valid bases, times 100
    GC_percentage = ((g_count + c_count) / valid_bases) * 100
    return GC_percentage


# ============================================
# STEP 4: CREATE A SIMPLE TABLE SHOWING GC CONTENT
# ============================================
# This just prints a quick table to the screen so we can see results immediately

def create_gc_summary_table(sequences_dict):
    """
    Makes a simple text table showing GC% for each strain.
    Think of it like a quick preview before we do the detailed report.
    """
    
    print("\n" + "="*70)
    print("GC CONTENT ANALYSIS SUMMARY TABLE")
    print("="*70)
    print(f"{'Strain':<20} {'GC Content (%)':>15} {'Genome Length':>15} {'N Count':>12} {'Valid Bases':>15}")
    print("-"*70)
    
    # Loop through each strain and calculate its GC%
    for name, seq in sequences_dict.items():
        gc = GC_content_calculator(seq)
        length = len(seq)
        n_count = seq.upper().count('N')  # Count how many unknown bases (N's)
        valid_bases = length - n_count    # Total minus unknowns = good bases
        
        print(f"{name:<20} {gc:>15.3f} {length:>15} {n_count:>12} {valid_bases:>15}")
    
    print("="*70)
    
    # Also show some statistics about all the strains together
    gc_values = [GC_content_calculator(seq) for seq in sequences_dict.values()]
    print(f"\nSUMMARY STATISTICS:")
    print(f"  Mean GC%: {sum(gc_values)/len(gc_values):.3f}%")
    print(f"  Range: {min(gc_values):.3f}% - {max(gc_values):.3f}%")
    print(f"  Standard Deviation: {np.std(gc_values):.4f}%")


# ============================================
# STEP 5: CREATE A DETAILED GC TABLE USING PANDAS
# ============================================
# Pandas makes nice tables that we can save as CSV files (Excel-like)

def create_pandas_gc_table(sequences_dict):
    """
    This creates a proper data table (DataFrame) with ALL the information.
    We'll use this later to save to CSV and make charts.
    """
    
    # Start with an empty list to collect data
    data = []
    
    # Loop through each strain
    for name, seq in sequences_dict.items():
        # Calculate GC content
        gc = GC_content_calculator(seq)
        
        # Get basic info about the sequence
        length = len(seq)
        seq_upper = seq.upper()
        
        # Count each base type
        a_count = seq_upper.count('A')
        t_count = seq_upper.count('T')
        g_count = seq_upper.count('G')
        c_count = seq_upper.count('C')
        n_count = seq_upper.count('N')
        
        # Calculate valid bases (total minus unknowns)
        valid_bases = length - n_count
        
        # Add all this info as a new row in our data list
        data.append({
            'Strain': name,
            'GC_Content_%': round(gc, 3),
            'Genome_Length': length,
            'A_Count': a_count,
            'T_Count': t_count,
            'G_Count': g_count,
            'C_Count': c_count,
            'N_Count': n_count,
            'Valid_Bases': valid_bases
        })
    
    # Turn our list into a proper pandas DataFrame (like an Excel spreadsheet)
    df = pd.DataFrame(data)
    
    # Get reference GC% (first row) to calculate differences
    ref_gc = df.iloc[0]['GC_Content_%']
    
    # Add a column showing how much each variant differs from reference
    df['Diff_from_Ref_%'] = df['GC_Content_%'] - ref_gc
    
    return df


# ============================================
# STEP 6: CREATE A FANCY, DETAILED GC REPORT
# ============================================
# This is the main GC analysis function - it does everything:
# 1. Makes a nice table
# 2. Calculates statistics
# 3. Saves to CSV
# 4. Ranks variants

def comprehensive_gc_report(sequences_dict):
    """
    This is our "all-in-one" GC analysis function.
    It prints a professional-looking report and saves everything to a file.
    """
    
    # First, get our data table using the function we just made
    df = create_pandas_gc_table(sequences_dict)
    
    # Print a nice header
    print("\n" + "="*90)
    print("COMPREHENSIVE GC CONTENT ANALYSIS REPORT")
    print("="*90)
    
    # Show the main results table
    print("\n📊 MAIN RESULTS:")
    
    # Pick only the columns we want to show
    display_columns = df[['Strain', 'GC_Content_%', 'Diff_from_Ref_%', 'Genome_Length', 'N_Count']]
    
    # Use tabulate to make it look fancy (like a formatted table)
    print(tabulate(display_columns, 
                   headers=['Strain', 'GC%', 'Δ from Ref', 'Length', 'Ns'],
                   tablefmt='fancy_grid', 
                   showindex=False,
                   floatfmt='.3f'))
    
    # Calculate and show summary statistics
    print("\n📈 SUMMARY STATISTICS:")
    
    stats = {
        'Mean GC%': df['GC_Content_%'].mean(),
        'Median GC%': df['GC_Content_%'].median(),
        'Std Dev': df['GC_Content_%'].std(),
        'Min GC%': df['GC_Content_%'].min(),
        'Max GC%': df['GC_Content_%'].max(),
        'Range': df['GC_Content_%'].max() - df['GC_Content_%'].min()
    }
    
    for stat_name, stat_value in stats.items():
        print(f"  {stat_name:15}: {stat_value:.4f}%")
    
    # Rank variants from highest to lowest GC%
    print("\n🏆 VARIANT RANKING (Highest to Lowest GC%):")
    
    # Sort by GC% (highest first)
    ranked = df.sort_values('GC_Content_%', ascending=False)
    
    # Give medals to top 3
    for i, row in ranked.iterrows():
        if i == 0:
            medal = "🥇"
        elif i == 1:
            medal = "🥈"
        elif i == 2:
            medal = "🥉"
        else:
            medal = "  "
        
        print(f"  {medal} {row['Strain']:25}: {row['GC_Content_%']:.3f}%")
    
    # Save everything to a CSV file
    try:
        # Try to save to home directory first
        home_path = os.path.expanduser("~/Complete_GC_Report.csv")
        df.to_csv(home_path, index=False)
        print(f"\n✅ Complete report saved to: {home_path}")
    except:
        # If that doesn't work, save to current folder
        df.to_csv('Complete_GC_Report.csv', index=False)
        print("\n✅ Complete report saved to: Complete_GC_Report.csv")
    
    return df


# ============================================
# STEP 7: FIND MUTATIONS IN EACH VARIANT
# ============================================
# This is the main mutation analysis function
# It compares each variant to the reference and finds all the differences

def analyze_all_variants(variants_dict, reference_seq, output_dir="mutation_results"):
    """
    Compare each variant to the reference sequence and find all mutations.
    Saves results to CSV files so we can analyze them later.
    """
    
    # Create a folder to store all mutation results (if it doesn't exist)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"  📁 Created folder: {output_dir}")
    
    all_mutations_data = {}  # Store all mutations for later use
    mutation_summary = []     # Store summary info (just counts)
    
    # Loop through each variant
    for variant_name, variant_seq in variants_dict.items():
        print(f"\n🔬 Analyzing {variant_name}...")
        
        # Set up the aligner (like putting two sequences side by side)
        aligner = PairwiseAligner()
        aligner.mode = "global"  # Compare the whole sequences, not just parts
        
        # Do the alignment
        alignments = aligner.align(reference_seq, variant_seq)
        alignment = alignments[0]  # Take the best alignment
        
        # Get the aligned sequences
        ref_aligned = alignment[0]  # Reference sequence after alignment
        var_aligned = alignment[1]  # Variant sequence after alignment
        
        # Find all positions where they differ
        mutations = []
        for i, (r, q) in enumerate(zip(ref_aligned, var_aligned)):
            # Skip if either base is unknown (N)
            if r == 'N' or q == 'N':
                continue
            # If they're different, that's a mutation!
            if r != q:
                mutations.append({
                    'Position': i+1,  # +1 because Python counts from 0, but biologists count from 1
                    'Reference': r,
                    'Variant': q,
                    'Variant_Strain': variant_name
                })
        
        # Save mutations to CSV file
        if mutations:
            df_mutations = pd.DataFrame(mutations)
            # Clean the filename (remove spaces)
            clean_name = variant_name.replace(' ', '_')
            filename = f"{output_dir}/{clean_name}_mutations.csv"
            df_mutations.to_csv(filename, index=False)
            print(f"  ✅ Found {len(mutations)} mutations - saved to {filename}")
            
            all_mutations_data[variant_name] = df_mutations
            mutation_summary.append({
                'Strain': variant_name,
                'Mutation_Count': len(mutations),
                # Mutation rate = mutations per 100 bases
                'Mutation_Rate': (len(mutations) / len(ref_aligned)) * 100
            })
        else:
            print(f"  ⚠️  No mutations found for {variant_name}")
            all_mutations_data[variant_name] = pd.DataFrame()
            mutation_summary.append({
                'Strain': variant_name,
                'Mutation_Count': 0,
                'Mutation_Rate': 0
            })
    
    # Save the summary (just counts) to a CSV file
    summary_df = pd.DataFrame(mutation_summary)
    summary_df.to_csv(f"{output_dir}/mutation_summary.csv", index=False)
    print(f"\n✅ Mutation summary saved to {output_dir}/mutation_summary.csv")
    
    return all_mutations_data, summary_df


# ============================================
# STEP 8: LOAD GENE INFORMATION FROM GENBANK FILE
# ============================================
# This reads the .gb file to find out where each gene is located

def load_gene_annotations(genbank_file="sequence.gb"):
    """
    Read the GenBank file to get information about gene locations.
    This tells us which part of the genome codes for which protein.
    """
    
    try:
        # Read the GenBank file
        record = SeqIO.read(genbank_file, "genbank")
        genes = []
        
        # Look through all the features in the genome
        for feature in record.features:
            # We only care about CDS = "Coding Sequence" (genes that make proteins)
            if feature.type == "CDS":
                # Get gene name (if available)
                gene_name = feature.qualifiers.get("gene", ["Unknown"])[0]
                # Get protein name (if available)
                product = feature.qualifiers.get("product", ["Unknown"])[0]
                # Get start and end positions
                start = int(feature.location.start)
                end = int(feature.location.end)
                
                genes.append({
                    'Gene': gene_name,
                    'Product': product,
                    'Start': start,
                    'End': end,
                    'Length': end - start
                })
        
        # Save to CSV for later use
        genes_df = pd.DataFrame(genes)
        genes_df.to_csv("gene_annotations.csv", index=False)
        print(f"✅ Loaded {len(genes)} genes from {genbank_file}")
        
        return genes_df
    except FileNotFoundError:
        print(f"⚠️  Gene annotation file '{genbank_file}' not found.")
        print("   Continuing without gene annotations...")
        return pd.DataFrame()


# ============================================
# STEP 9: FIGURE OUT WHICH GENE EACH MUTATION IS IN
# ============================================
# This takes each mutation and checks what gene it falls in

def annotate_mutations_with_genes(mutations_dict, genes_df):
    """
    For each mutation, figure out if it's inside a gene or between genes.
    If it's inside a gene, record which gene.
    """
    
    if genes_df.empty:
        print("⚠️  No gene annotations available. Skipping annotation.")
        return mutations_dict
    
    annotated_mutations = {}
    
    # Loop through each variant's mutations
    for strain, mutations_df in mutations_dict.items():
        if mutations_df.empty:
            annotated_mutations[strain] = mutations_df
            continue
        
        # Make a copy so we don't mess up the original
        df_annotated = mutations_df.copy()
        
        # Add two new columns
        df_annotated['Gene'] = 'Intergenic'  # Default: between genes
        df_annotated['Product'] = 'None'      # Default: no protein
        
        # Check each mutation position
        for idx, row in df_annotated.iterrows():
            pos = row['Position']
            
            # Look through all genes to see if this position falls inside one
            for _, gene in genes_df.iterrows():
                if gene['Start'] <= pos <= gene['End']:
                    df_annotated.at[idx, 'Gene'] = gene['Gene']
                    df_annotated.at[idx, 'Product'] = gene['Product']
                    break  # Stop looking once we find the gene
        
        # Save annotated mutations
        clean_name = strain.replace(' ', '_')
        filename = f"mutation_results/{clean_name}_annotated_mutations.csv"
        df_annotated.to_csv(filename, index=False)
        annotated_mutations[strain] = df_annotated
    
    return annotated_mutations


# ============================================
# STEP 10: CREATE BAR CHART OF GC CONTENT
# ============================================
# This makes Figure 1 - a bar chart comparing GC% across variants

def plot_gc_content(sequences_dict, output_file="gc_content_bar_chart.png"):
    """
    Make a bar chart showing GC content for each strain.
    This helps visualize which variants have higher/lower GC%.
    """
    
    strains = []
    gc_values = []
    colors = ['#2E86AB', '#A23B72', '#F18F01', '#C73E1D', '#6A994E']  # Nice colors
    
    # Get GC% for each strain
    for name, seq in sequences_dict.items():
        strains.append(name)
        gc_values.append(GC_content_calculator(seq))
    
    # Create the figure
    plt.figure(figsize=(12, 6))
    bars = plt.bar(strains, gc_values, color=colors, edgecolor='black', linewidth=1)
    
    # Add labels and title
    plt.ylabel('GC Content (%)', fontsize=12)
    plt.title('Figure 1: GC Content Comparison Across SARS-CoV-2 Variants', 
              fontsize=14, fontweight='bold')
    plt.ylim(37.85, 38.05)  # Zoom in to see small differences
    plt.xticks(rotation=45, ha='right')  # Rotate strain names for readability
    plt.grid(axis='y', alpha=0.3, linestyle='--')
    
    # Add the actual numbers on top of each bar
    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height + 0.01,
                f'{height:.3f}%', ha='center', va='bottom', fontweight='bold')
    
    # Add a reference line showing the reference strain's GC%
    plt.axhline(y=gc_values[0], color='red', linestyle='--', 
                linewidth=1.5, alpha=0.7, label=f"Reference: {gc_values[0]:.3f}%")
    plt.legend()
    
    # Save and display
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.show()
    print(f"✅ GC content chart saved as {output_file}")
    
    # Save the figure caption for our report
    with open("figure_captions.txt", "a") as f:
        f.write("\nFigure 1: Bar chart comparing GC content percentages across five SARS-CoV-2 "
                "variants (Wuhan-Hu-1 reference, Alpha, Delta, Omicron BA.1, and Omicron XBB). "
                f"Values range from {min(gc_values):.3f}% to {max(gc_values):.3f}%.\n")


# ============================================
# STEP 11: CREATE BAR CHART OF MUTATION COUNTS
# ============================================
# This makes Figure 2 - showing how many mutations each variant has

def plot_mutation_counts(mutation_summary_df, output_file="mutation_count_chart.png"):
    """
    Make a bar chart showing the total number of mutations in each variant.
    Higher bars = more mutations = more evolved/different from original.
    """
    
    plt.figure(figsize=(10, 6))
    
    strains = mutation_summary_df['Strain']
    counts = mutation_summary_df['Mutation_Count']
    
    # Create bars with different colors
    bars = plt.bar(strains, counts, color=['#FF6B6B', '#4ECDC4', '#45B7D1', '#96CEB4'])
    
    # Labels and title
    plt.ylabel('Number of Mutations', fontsize=12)
    plt.xlabel('Variant Strain', fontsize=12)
    plt.title('Figure 2: Mutation Counts Relative to Wuhan-Hu-1 Reference', 
              fontsize=14, fontweight='bold')
    plt.xticks(rotation=45, ha='right')
    plt.grid(axis='y', alpha=0.3, linestyle='--')
    
    # Add numbers on top of bars
    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height + 0.5,
                f'{int(height)}', ha='center', va='bottom', fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.show()
    print(f"✅ Mutation count chart saved as {output_file}")
    
    # Save caption
    with open("figure_captions.txt", "a") as f:
        f.write(f"\nFigure 2: Bar chart showing the total number of point mutations identified "
                f"in each variant compared to the Wuhan-Hu-1 reference genome.\n")


# ============================================
# STEP 12: CREATE MUTATION DISTRIBUTION PLOT
# ============================================
# This makes Figure 3 - shows WHERE mutations occur in the genome
# Each dot is a mutation, so we can see which regions have lots of mutations

def plot_mutation_distribution(annotated_mutations_dict, genes_df, output_file="mutation_distribution.png"):
    """
    Make a scatter plot showing the position of each mutation along the genome.
    This helps us see mutation hotspots visually.
    """
    
    plt.figure(figsize=(15, 8))
    
    colors = ['red', 'blue', 'green', 'purple']
    y_pos = 0  # Starting position for first variant
    
    # Plot mutations for each strain at different vertical positions
    for i, (strain, mutations_df) in enumerate(annotated_mutations_dict.items()):
        if mutations_df.empty:
            continue
            
        positions = mutations_df['Position'].values
        # Plot each mutation as a dot, slightly offset vertically so they don't overlap
        plt.scatter(positions, [y_pos + i*0.1] * len(positions), 
                   c=colors[i % len(colors)], s=30, alpha=0.6, 
                   label=strain, edgecolors='black', linewidth=0.5)
    
    # Add gray boxes showing where genes are located
    if not genes_df.empty:
        for _, gene in genes_df.iterrows():
            plt.axvspan(gene['Start'], gene['End'], alpha=0.1, color='gray')
            plt.text(gene['Start'] + (gene['Length']/2), -0.2, 
                    gene['Gene'], rotation=45, fontsize=8, ha='center')
    
    plt.xlabel('Genome Position (bp)', fontsize=12)
    plt.ylabel('Strains', fontsize=12)
    plt.title('Figure 3: Distribution of Mutations Across the SARS-CoV-2 Genome', 
              fontsize=14, fontweight='bold')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True, alpha=0.3, linestyle='--')
    plt.tight_layout()
    
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.show()
    print(f"✅ Mutation distribution plot saved as {output_file}")
    
    with open("figure_captions.txt", "a") as f:
        f.write(f"\nFigure 3: Scatter plot showing the distribution of mutations across the "
                f"SARS-CoV-2 genome for each variant. Gray shaded regions indicate gene coding "
                f"sequences.\n")


# ============================================
# STEP 13: CREATE MUTATION HEATMAP
# ============================================
# This makes Figure 4 - a heatmap showing mutation density
# Darker colors = more mutations in that region

def create_mutation_heatmap(annotated_mutations_dict, genome_length=29903, 
                           bin_size=500, output_file="mutation_heatmap.png"):
    """
    Make a heatmap showing mutation density across the genome.
    We divide the genome into 500-base chunks and count mutations in each chunk.
    """
    
    # Create bins (chunks) across the genome
    bins = np.arange(0, genome_length + bin_size, bin_size)
    
    # Calculate mutation density for each strain in each bin
    heatmap_data = []
    strain_names = []
    
    for strain, mutations_df in annotated_mutations_dict.items():
        strain_names.append(strain)
        if mutations_df.empty:
            # If no mutations, all bins have 0
            density = [0] * (len(bins) - 1)
        else:
            positions = mutations_df['Position'].values
            # Count how many mutations fall into each bin
            density, _ = np.histogram(positions, bins=bins)
        heatmap_data.append(density)
    
    # Create a table (DataFrame) for the heatmap
    heatmap_df = pd.DataFrame(heatmap_data, 
                             index=strain_names,
                             columns=[f"{int(bins[i])}-{int(bins[i+1])}" for i in range(len(bins)-1)])
    
    # Draw the heatmap
    plt.figure(figsize=(16, 6))
    sns.heatmap(heatmap_df, cmap='YlOrRd', annot=False, cbar_kws={'label': 'Mutation Count'})
    # YlOrRd = Yellow to Orange to Red - darker red = more mutations
    
    plt.xlabel('Genome Position (bp)', fontsize=12)
    plt.ylabel('Variant Strain', fontsize=12)
    plt.title('Figure 4: Mutation Density Heatmap Across SARS-CoV-2 Genome', 
              fontsize=14, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.show()
    print(f"✅ Mutation heatmap saved as {output_file}")
    
    with open("figure_captions.txt", "a") as f:
        f.write(f"\nFigure 4: Heatmap showing mutation density across the SARS-CoV-2 genome "
                f"(binned in 500bp intervals). Darker red colors indicate higher mutation "
                f"frequency.\n")


# ============================================
# STEP 14: ORGANIZE ALL OUR SEQUENCES
# ============================================
# Put all sequences into dictionaries so we can loop through them easily

# All sequences (including reference) for GC analysis
sequences = {
    "Reference (Wuhan-Hu-1)": reference_sequence,
    "Alpha Variant": alpha_variant,
    "Delta Variant": delta_variant,
    "Omicron BA.1": Omicron_BA_1,
    "Omicron XBB": Omicron_XBB
}

# Only variants (no reference) for mutation analysis
variants_only = {
    "Alpha Variant": alpha_variant,
    "Delta Variant": delta_variant,
    "Omicron BA.1": Omicron_BA_1,
    "Omicron XBB": Omicron_XBB
}


# ============================================
# STEP 15: RUN THE ENTIRE ANALYSIS
# ============================================
# This is where everything actually happens!
# We call all the functions we defined above

print("\n" + "="*80)
print("🦠 SARS-CoV-2 COMPREHENSIVE GENOMIC ANALYSIS PIPELINE")
print("="*80)

# ===== PHASE 1: GC CONTENT ANALYSIS =====
print("\n📊 PHASE 1: GC CONTENT ANALYSIS")
print("-"*50)

# Show quick summary table
create_gc_summary_table(sequences)

# Show detailed report and save to CSV
results_df = comprehensive_gc_report(sequences)

# Create Figure 1 - GC content bar chart
plot_gc_content(sequences)


# ===== PHASE 2: MUTATION ANALYSIS =====
print("\n🧬 PHASE 2: MUTATION ANALYSIS")
print("-"*50)

# Find all mutations in each variant
mutations_dict, mutation_summary = analyze_all_variants(
    variants_only, 
    reference_sequence,
    output_dir="mutation_results"
)

# Create Figure 2 - Mutation count bar chart
if not mutation_summary.empty:
    plot_mutation_counts(mutation_summary)


# ===== PHASE 3: GENE ANNOTATION =====
print("\n📝 PHASE 3: GENE ANNOTATION")
print("-"*50)

# Load gene locations from GenBank file
genes_df = load_gene_annotations("sequence.gb")

# If we have gene info, do the advanced analysis
if not genes_df.empty:
    # Figure out which gene each mutation is in
    annotated_mutations = annotate_mutations_with_genes(mutations_dict, genes_df)
    
    # Create Figure 3 - Mutation distribution plot
    plot_mutation_distribution(annotated_mutations, genes_df)
    
    # Create Figure 4 - Mutation heatmap
    create_mutation_heatmap(annotated_mutations)
else:
    print("⚠️  Skipping gene-based visualizations (no gene annotation file found).")


# ===== FINAL SUMMARY =====
print("\n" + "="*80)
print("✅ ANALYSIS COMPLETE! Here's what was created:")
print("="*80)
print("📊 CSV Files:")
print("  - Complete_GC_Report.csv (GC content summary)")
print("  - mutation_results/mutation_summary.csv (Mutation counts)")
print("  - mutation_results/*_mutations.csv (Individual variant mutations)")
if os.path.exists("gene_annotations.csv"):
    print("  - gene_annotations.csv (Gene positions)")
print("\n🖼️  Figures:")
print("  - gc_content_bar_chart.png (Figure 1)")
print("  - mutation_count_chart.png (Figure 2)")
if os.path.exists("gene_annotations.csv"):
    print("  - mutation_distribution.png (Figure 3)")
    print("  - mutation_heatmap.png (Figure 4)")
print("\n📝 Captions:")
print("  - figure_captions.txt (All figure captions)")
print("="*80)