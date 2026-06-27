#!/usr/bin/env python3
import argparse
import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings('ignore', r'All-NaN (slice|axis) encountered')
warnings.filterwarnings('ignore', r'Mean of empty slice')
import plotly.express as px
import plotly.graph_objects as go
import os

def parse_args():
    parser = argparse.ArgumentParser(description="Analyze and compare alignment strategies.")
    parser.add_argument("--events-tsv", required=True, help="Annotated events TSV file")
    parser.add_argument("--out-html", required=True, help="Path to output HTML report")
    return parser.parse_args()

def calc_titv(df):
    snvs = df[df['event_type'] == 'snv']
    if len(snvs) == 0:
        return np.nan
    
    transitions = {('A', 'G'), ('G', 'A'), ('C', 'T'), ('T', 'C')}
    
    ti = 0
    tv = 0
    for _, row in snvs.iterrows():
        ref, alt = row['ref'], row['alt']
        if len(ref) == 1 and len(alt) == 1:
            if (ref, alt) in transitions:
                ti += 1
            else:
                tv += 1
    
    if tv == 0:
        return np.nan if ti == 0 else float('inf')
    return round(ti / tv, 3)

def categorize_clinvar(val):
    if pd.isna(val) or val == "": return "Not Found"
    val = str(val).lower()
    if "pathogenic" in val: return "P/LP"
    if "benign" in val: return "B/LB"
    if "uncertain" in val: return "VUS"
    return "Other"

def main():
    args = parse_args()
    
    print(f"Reading {args.events_tsv}...")
    df = pd.read_csv(args.events_tsv, sep="\t", compression="gzip", low_memory=False)
    
    print("Preprocessing data...")
    df['clinvar_category'] = df['clinvar_sig'].apply(categorize_clinvar)
    df['gnomad_af'] = pd.to_numeric(df['gnomad_af'], errors='coerce')
    
    # Define a human variant ID
    df['variant_id'] = df['genomic_accession'].astype(str) + ":" + df['genomic_start1'].astype(str) + ":" + df['ref'].fillna('') + ">" + df['alt'].fillna('')
    
    strategies = df['strategy'].unique()
    html_sections = []
    
    print("Computing global metrics...")
    global_stats = []
    strategy_variant_sets = {}
    
    for strategy in strategies:
        s_df = df[df['strategy'] == strategy]
        unique_vars = s_df.drop_duplicates('variant_id').copy()
        strategy_variant_sets[strategy] = set(unique_vars['variant_id'])
        
        support_counts = s_df.groupby('variant_id')['ortholog_gene_id'].nunique()
        unique_vars['ortholog_count'] = unique_vars['variant_id'].map(support_counts)
        
        ti_tv = calc_titv(unique_vars)
        clinvar_counts = unique_vars['clinvar_category'].value_counts()
        
        global_stats.append({
            'Strategy': strategy,
            'Total Unique Variants': len(unique_vars),
            'Ti/Tv Ratio': ti_tv,
            'Pathogenic (P/LP)': clinvar_counts.get('P/LP', 0),
            'Benign (B/LB)': clinvar_counts.get('B/LB', 0),
            'VUS': clinvar_counts.get('VUS', 0),
            'Median gnomAD AF': unique_vars['gnomad_af'].median(),
            'Mean gnomAD AF': unique_vars['gnomad_af'].mean()
        })
        
    global_df = pd.DataFrame(global_stats)
    html_sections.append("<h2>Global Strategy Comparison</h2>")
    html_sections.append(global_df.to_html(index=False, classes='table table-striped table-bordered', float_format='%.5g'))
    
    print("Computing ortholog support breakdown...")
    html_sections.append("<h2>Metrics by Ortholog Support Count</h2>")
    html_sections.append("<p>Variants grouped by how many distinct orthologs aligned to produce them.</p>")
    
    for strategy in strategies:
        s_df = df[df['strategy'] == strategy]
        unique_vars = s_df.drop_duplicates('variant_id').copy()
        support_counts = s_df.groupby('variant_id')['ortholog_gene_id'].nunique()
        unique_vars['ortholog_count'] = unique_vars['variant_id'].map(support_counts)
        
        unique_vars['ortholog_bucket'] = unique_vars['ortholog_count'].apply(lambda x: str(x) if x < 5 else "5+")
        
        stats = []
        buckets = sorted(unique_vars['ortholog_bucket'].unique(), key=lambda x: int(x.replace('+','')))
        for b in buckets:
            b_df = unique_vars[unique_vars['ortholog_bucket'] == b]
            clin_c = b_df['clinvar_category'].value_counts()
            stats.append({
                'Ortholog Support': b,
                'Variant Count': len(b_df),
                'Ti/Tv': calc_titv(b_df),
                'P/LP': clin_c.get('P/LP', 0),
                'B/LB': clin_c.get('B/LB', 0),
                'VUS': clin_c.get('VUS', 0),
                'Median AF': b_df['gnomad_af'].median()
            })
            
        b_df_stats = pd.DataFrame(stats)
        html_sections.append(f"<h3>{strategy}</h3>")
        html_sections.append(b_df_stats.to_html(index=False, classes='table table-sm', float_format='%.5g'))
        
    print("Computing unique contributions...")
    html_sections.append("<h2>Unique Variants per Strategy</h2>")
    html_sections.append("<p>Variants found ONLY by this strategy.</p>")
    
    unique_stats = []
    for strategy in strategies:
        other_strats = [s for s in strategies if s != strategy]
        union_others = set()
        for os_s in other_strats:
            union_others.update(strategy_variant_sets[os_s])
            
        my_unique = strategy_variant_sets[strategy] - union_others
        
        if len(my_unique) > 0:
            s_df = df[(df['strategy'] == strategy) & (df['variant_id'].isin(my_unique))].drop_duplicates('variant_id')
            clin_c = s_df['clinvar_category'].value_counts()
            
            unique_stats.append({
                'Strategy': strategy,
                'Unique Variants': len(my_unique),
                'Ti/Tv': calc_titv(s_df),
                'P/LP': clin_c.get('P/LP', 0),
                'B/LB': clin_c.get('B/LB', 0),
                'VUS': clin_c.get('VUS', 0),
                'Median AF': s_df['gnomad_af'].median()
            })
        else:
            unique_stats.append({
                'Strategy': strategy,
                'Unique Variants': 0
            })
            
    html_sections.append(pd.DataFrame(unique_stats).to_html(index=False, classes='table table-bordered', float_format='%.5g'))
    
    print("Generating plots...")
    html_sections.append("<h2>Distributions</h2>")
    
    plot_df = df.drop_duplicates(subset=['strategy', 'variant_id'])
    
    fig1 = px.histogram(plot_df[plot_df['clinvar_category'] != 'Not Found'], 
                        x='strategy', color='clinvar_category', barmode='group',
                        title="ClinVar Classifications Found by Strategy (excluding Not Found)",
                        category_orders={"clinvar_category": ["P/LP", "B/LB", "VUS", "Other"]})
    html_sections.append(fig1.to_html(full_html=False, include_plotlyjs='cdn'))
    
    # AF Boxplot
    # use log10 AF for better visibility
    plot_df = plot_df.copy()
    plot_df['log10_gnomad_af'] = np.log10(plot_df['gnomad_af'].replace(0, np.nan))
    fig2 = px.box(plot_df.dropna(subset=['log10_gnomad_af']), 
                  x='strategy', y='log10_gnomad_af', 
                  title="gnomAD AF Distribution (Log10) by Strategy",
                  points="outliers")
    html_sections.append(fig2.to_html(full_html=False, include_plotlyjs='cdn'))

    print(f"Writing report to {args.out_html}...")
    html_template = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <meta charset="utf-8">
        <title>Strategy Comparison Report</title>
        <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/css/bootstrap.min.css" rel="stylesheet">
        <style>
            body {{ padding: 20px; font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif; }}
            h2 {{ margin-top: 40px; border-bottom: 1px solid #ccc; padding-bottom: 10px; }}
            h3 {{ margin-top: 20px; }}
            .table {{ width: auto; max-width: 100%; margin-bottom: 1rem; background-color: transparent; }}
        </style>
    </head>
    <body>
        <div class="container-fluid">
            <h1>Alignment Strategies Report</h1>
            <p class="lead">Comprehensive comparison of variants discovered by different alignment strategies.</p>
            {''.join(html_sections)}
        </div>
    </body>
    </html>
    """
    
    with open(args.out_html, "w") as f:
        f.write(html_template)
        
    print("Done!")

if __name__ == "__main__":
    main()
