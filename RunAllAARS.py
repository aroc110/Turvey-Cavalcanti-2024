import pandas as pd
from datetime import datetime
import numpy as np
import plotly.express as px
from plotly import graph_objects as go

COVERAGE = 5 # Minimum coverage of the alignment to be considered present in bacteria or archaea
MINDOMAINLEN = 15 # Minimum number of residues in a domain to be considered modern
## Constant parameters
# LINEAGES -- dictionary with the taxonomic lineage IDs for the taxonomic groups to be included in the analysis.
# ORDERTAX -- list with the taxonomic groups to be included in the analysis in order of least specific to most specific when one group is contained in another. Here vertebrates need to come before mammals.
LINEAGES = {
    "Bacteria": 2,
    "Archaea": 2157,
    "Vertebrata": 7742,
    "Mammals": 40674
}
ORDERTAX = ["Bacteria", "Archaea", "Vertebrata", "Mammals"]
WANTEDTAX = ["Mammals","Vertebrata"]

ECNUMBERS = {
    "Ala" : ['6.1.1.7', 'SYAC_HUMAN', 'AARS1'],
    "Arg" : ['6.1.1.19', 'SYRC_HUMAN', 'RARS1'],
    "Asn" : ['6.1.1.22', 'SYNC_HUMAN', 'NARS1'],
    "Asp" : ['6.1.1.12', 'SYDC_HUMAN', 'DARS1'],
    "Cys" : ['6.1.1.16', 'SYCC_HUMAN', 'CARS1'],

    "Gln" : ['6.1.1.18', 'SYQ_HUMAN', 'QARS1'],
    "Gly" : ['6.1.1.14', 'GARS_HUMAN', 'GARS1'],
    "His" : ['6.1.1.21', 'HARS1_HUMAN', 'HARS1'],
    "Ile" : ['6.1.1.5', 'SYIC_HUMAN', 'IARS1'],

    "Leu" : ['6.1.1.4', 'SYLC_HUMAN', 'LARS1'],
    "Lys" : ['6.1.1.6', 'SYK_HUMAN', 'KARS1'],
    "Met" : ['6.1.1.10', 'SYMC_HUMAN', 'MARS1'],

    'PheA' : ['6.1.1.20', 'SYFA_HUMAN', 'FARSA'],
    'PheB' : ['6.1.1.20', 'SYFB_HUMAN', 'FARSB'], # This is a special case where two genes have the same EC number they are analyzed separately and we combine results manually in end.
    "Ser" : ['6.1.1.11', 'SYSC_HUMAN', 'SARS1'],

    "Thr" : ['6.1.1.3', 'SYTC_HUMAN', 'TARS1'],
    "Trp" : ['6.1.1.2', 'SYWC_HUMAN', 'WARS1'],
    "Tyr" : ['6.1.1.1', 'SYYC_HUMAN', 'YARS1'],
    "Val" : ['6.1.1.9', 'SYVC_HUMAN', 'VARS1'],
    "GluPro" : [['6.1.1.17','6.1.1.15'], 'SYEP_HUMAN', 'EPRS1'], # The ProRS and GluRS genes are fused
}

def write_report(name, coverage, agreement, uTestsCons, conservation1, conservation2, stats, tables, tableDomains):
    fileout = open(f"{name}/{name}.html","w")
    fileout.write(f'''<!DOCTYPE html>
<html>
<head>
    <title>Home</title>
    <style>
        .grid-1 {{
          display: grid;
          grid-template-columns: 1fr 1fr;
          grid-template-rows: 1fr;
        }}
        .table2 {{ grid-area: 1 / 1 / 2 / 2; }}
        .table3 {{ grid-area: 1 / 2 / 2 / 3; }}
        
        .grid-container {{
          display: grid;
          grid-template-columns: 1fr 1.5fr;
          grid-template-rows: 1fr 1fr 1.5fr;
        }}
        h1, h2, h3, p {{
            font-family: Arial, sans-serif;
            }}
        p {{
            width: 80%;
            margin: auto;
        }}
        .coverage {{ grid-area: 1 / 1 / 2 / 3; }}
        .confusion {{ grid-area: 2 / 1 / 3 / 2; }}
        .significance {{ grid-area: 2 / 2 / 3 / 3; }}
        .conservation {{ grid-area: 3 / 1 / 4 / 3; }}
        .pandas-table {{
            font-family: Arial, Helvetica, sans-serif;
            border-collapse: collapse;
            width: 80%;
            margin: auto;
        }}

        .pandas-table td, .pandas-table th {{
            border: 1px solid #ddd;
            padding: 8px;
            text-align: center;
        }}

        .pandas-table tr:nth-child(even){{
            background-color: #f2f2f2;}}

        .pandas-table tr:hover {{
            background-color: #ddd;
        }}

        .pandas-table th {{
            padding-top: 12px;
            padding-bottom: 12px;
            text-align: left;
            background-color: #0351ab;
            color: white;
        }}

        </style>
</head>
<body>
    <h1>Analysis of {name} aminoacyl synthetase.</h1>
    <p>Analysis performed on {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}</p>
    <h2>Number of Proteins Analyzed</h2>
    {stats.to_html(classes="pandas-table")}
    
    <h2>Number of Variants Analyzed:</h2>
    {tables[0].to_html(classes="pandas-table",index=False)}
    
    <h2>Domains</h2>
    {tableDomains.to_html(classes="pandas-table",float_format="{:.0f}".format)}
    
    <h2>Results</h2>
    <h3>Mutations in Ancient vs Modern Domains</h3>
    <div class="grid-1">
        <div class="table2">
        <h3>Current Study</h3>
        {tables[1].to_html(classes="pandas-table")}
        <p>The Fisher exact test was calculated using a contingency table with the number of pathogenic and benign/unknown variants in ancient and modern domains. The bottom left values.</p>
        </div>
        <div class="table3">
        <h3>Guo et al (2010)</h3>
        {tables[2].to_html(classes="pandas-table")}
        <p>The Fisher exact test was calculated using a contingency table with the number of pathogenic and benign/unknown variants in ancient and modern domains. The bottom left values.</p>
        </div>
    </div>

    <div class="grid-container">
        <div class="coverage">
            <h3>Coverage of the alignment and domain locations</h2>
            {coverage.to_html()}
        </div>
        <div class="confusion">
            <h3>Agreement between the current study and Guo et al (2010)</h2>
            {agreement.to_html()}
        </div>
        <div class="significance">
            <h3>Conservation of the amino acids in modern and ancient domains</h2>
            {uTestsCons.to_html(classes="pandas-table",float_format="{:.3e}".format)}
            <p>The Mann-Whitney U test was calculated using the conservation scores of the amino acids in ancient and modern domains. The p-value is for the null hypothesis that the conservation scores in ancient domains are greater than in modern domains.</p>
        </div>
        <div class="conservation">
            <h3>Conservation of the amino acids in modern and ancient domains</h2>
            {conservation1.to_html()}
            {conservation2.to_html()}
        </div>
      </div>

</body>
</html>''')
    fileout.close()
    return(0)

def make_test_table(resultGraph, variants):
    newDf = pd.merge(variants, resultGraph[['Position','Current Study','Guo et al (2010)']],  how='left', left_on=['Pos'], right_on = ['Position'])

    # Table 1 with the number of pathogenic and benign/unknown variants in the full protein
    protlength = len(resultGraph)
    protLen = len(resultGraph)
    pathogenic = len(variants[variants['Pathogenicity'] == 'Pathogenic'])
    benign = len(variants[variants['Pathogenicity'] == 'Benign/Unknown'])
    
    tlist1 = [[protLen, pathogenic, benign]]
    table1 = pd.DataFrame(tlist1, columns=["Protein Length", "Pathogenic Variants", "Benign/Unknown Variants"]) 
    
    # Table 2: Compare the number of variants in ancient vs modern domains in the current study

    fulltable = [["count","Ancient Domains","Modern Domains","All","Ratio"]]

    # Number of residues
    ancientAll = len(resultGraph[resultGraph['Current Study'] == 'Ancient'])
    modernAll = len(resultGraph[resultGraph['Current Study'] == 'Modern'])
    allAll = ancientAll + modernAll
    ratioAll = ancientAll/modernAll if modernAll > 0 else "Undefined"
    fulltable.append(["Number of Residues", ancientAll, modernAll, allAll, ratioAll])

    # Number of variants
    allMut = len(newDf)
    ancientMut = len(newDf[newDf['Current Study'] == 'Ancient'])
    modernMut = len(newDf[newDf['Current Study'] == 'Modern'])
    ratioMut = ancientMut/modernMut if modernMut > 0 else "Undefined"
    fulltable.append(["Number of Variants", ancientMut, modernMut, allMut, ratioMut])

    # Number of pathogenic variants
    allPatho = len(newDf[newDf['Pathogenicity'] == 'Pathogenic'])
    ancientPatho = len(newDf[(newDf['Current Study'] == 'Ancient') & (newDf['Pathogenicity'] == 'Pathogenic')])
    modernPatho = len(newDf[(newDf['Current Study'] == 'Modern') & (newDf['Pathogenicity'] == 'Pathogenic')])
    ratioPatho = ancientPatho/modernPatho if modernPatho > 0 else "Undefined"
    fulltable.append(["Number of Pathogenic Variants", ancientPatho, modernPatho, allPatho, ratioPatho])

    # Number of benign/unknown variants
    allUnk = len(newDf[newDf['Pathogenicity'] == 'Benign/Unknown'])
    ancientUnk = len(newDf[(newDf['Current Study'] == 'Ancient') & (newDf['Pathogenicity'] == 'Benign/Unknown')])
    modernUnk = len(newDf[(newDf['Current Study'] == 'Modern') & (newDf['Pathogenicity'] == 'Benign/Unknown')])
    ratioUnk = ancientUnk/modernUnk if modernUnk > 0 else "Undefined"
    fulltable.append(["Number of Benign/Unknown Variants", ancientUnk, modernUnk, allUnk, ratioUnk])

    contingency = np.array([
    [ancientPatho, modernPatho],
    [ancientUnk, modernUnk]
    ])

    from scipy.stats import fisher_exact as fe
    if modernPatho + modernUnk == 0:
        fisher = [0,"Undefined"]
    else:
        fisher = fe(contingency,alternative="two-sided")
    fulltable.append(["Fisher exact test","","","",fisher[1]])

    table2 = pd.DataFrame(fulltable[1:], columns=fulltable[0])
    table2.set_index("count", inplace=True)
    table2.index.name = None
    
    # Table 3: Compare the number of variants in ancient vs modern domains using Guo et al (2010) data

    fulltable3 = [["count","Ancient Domains","Modern Domains","All","Ratio"]]

    # Number of residues
    ancientAll3 = len(resultGraph[resultGraph['Guo et al (2010)'] == 'Ancient'])
    modernAll3 = len(resultGraph[resultGraph['Guo et al (2010)'] == 'Modern'])
    allAll3 = ancientAll3 + modernAll3
    ratioAll3 = ancientAll3/modernAll3 if modernAll3 > 0 else "Undefined"
    fulltable3.append(["Number of Residues", ancientAll3, modernAll3, allAll3, ratioAll3])

    # Number of variants
    allMut3 = len(newDf)
    ancientMut3 = len(newDf[newDf['Guo et al (2010)'] == 'Ancient'])
    modernMut3 = len(newDf[newDf['Guo et al (2010)'] == 'Modern'])
    ratioMut3 = ancientMut3/modernMut3 if modernMut3 > 0 else "Undefined"
    fulltable3.append(["Number of Variants", ancientMut3, modernMut3, allMut3, ratioMut3])

    # Number of pathogenic variants
    allPatho3 = len(newDf[newDf['Pathogenicity'] == 'Pathogenic'])
    ancientPatho3 = len(newDf[(newDf['Guo et al (2010)'] == 'Ancient') & (newDf['Pathogenicity'] == 'Pathogenic')])
    modernPatho3 = len(newDf[(newDf['Guo et al (2010)'] == 'Modern') & (newDf['Pathogenicity'] == 'Pathogenic')])
    ratioPatho3 = ancientPatho3/modernPatho3 if modernPatho3 > 0 else "Undefined"
    fulltable3.append(["Number of Pathogenic Variants", ancientPatho3, modernPatho3, allPatho3, ratioPatho3])

    # Number of benign/unknown variants
    allUnk3 = len(newDf[newDf['Pathogenicity'] == 'Benign/Unknown'])
    ancientUnk3 = len(newDf[(newDf['Guo et al (2010)'] == 'Ancient') & (newDf['Pathogenicity'] == 'Benign/Unknown')])
    modernUnk3 = len(newDf[(newDf['Guo et al (2010)'] == 'Modern') & (newDf['Pathogenicity'] == 'Benign/Unknown')])
    ratioUnk3 = ancientUnk3/modernUnk3 if modernUnk3 > 0 else "Undefined"
    fulltable3.append(["Number of Benign/Unknown Variants", ancientUnk3, modernUnk3, allUnk3, ratioUnk3])

    contingency3 = np.array([
    [ancientPatho3, modernPatho3],
    [ancientUnk3, modernUnk3]
    ])

    from scipy.stats import fisher_exact as fe
    if modernPatho3 + modernUnk3 == 0:
        fisher3 = [0,"Undefined"]
    else:
        fisher3 = fe(contingency3,alternative="two-sided")
    fulltable3.append(["Fisher exact test","","","",fisher3[1]])

    table3 = pd.DataFrame(fulltable3[1:], columns=fulltable3[0])
    table3.set_index("count", inplace=True)
    table3.index.name = None

    return([table1, table2, table3])

def make_coverage_plot(coverage, annotation,aa):
    '''Makes a plot showing the taxon coverage of the alignment and domain locations.
    
    Arguments:
    coverage -- a dataframe with the coverage information
    annotations -- a list of lists with the domain annotations
    
    '''
    from plotly.subplots import make_subplots

    rowHeights = [0.9,0.1]
    fig = make_subplots(rows=2, cols=1, shared_xaxes=True, row_heights=rowHeights, vertical_spacing=0.02)


    ys = ["Bacteria %", "Archaea %", "Vertebrata %", "Mammals %"]
    for yaxis in ys:
        fig.add_trace(go.Scatter(x=coverage['Position'], y=coverage[yaxis],name=yaxis, legendgroup='Taxon Coverage',legendgrouptitle_text="Taxon Coverage"),row=1,col=1)


    for i in annotation:
        fig.add_shape(
            type="rect",
            xref="x", yref="y",
            x0=i[0], y0=0, x1=i[1], y1=100,
            opacity=0.2,
            fillcolor="orange",
            line_color="orange",
        )
        

    ## Add modern domains
    newdoms = list(coverage['Current Study'].apply(lambda x: 0 if x == 'Ancient' else 1))
    
    fig.add_trace(go.Bar(x=coverage['Position'], y=newdoms, name='Current Study', marker_color='blue', marker_line_width=0, marker_line_color='blue', width=1, legendgroup='Modern Domains',legendgrouptitle_text="Modern Domains"), row=2, col=1)
    
    fig.add_trace(go.Bar(x=coverage['Position'], y=[0]*len(coverage['Position']), name='Guo et al (2010)', marker_color='orange', marker_line_width=0, marker_line_color='orange', width=0,opacity=0.2, legendgroup='Modern Domains'), row=2, col=1)



    fig.update_layout(template="simple_white",title=aa) #,width=1200,height=400)



    fig.update_yaxes(range=[0,100], row=1,col=1)
    fig.update_yaxes(tickvals=[0,100],ticks='',row=1,col=1)
    fig.update_xaxes(ticks='', row=1,col=1)
    fig.update_yaxes(zeroline=True, row=1,col=1)
    fig.update_xaxes(showline=False, row=1,col=1)


    fig.update_yaxes(ticks='', showticklabels=False, row=2,col=1)
    fig.update_xaxes(ticks='', row=2,col=1)



    fig.update_xaxes(title='Amino Acid Position', row=2,col=1)
    fig.update_yaxes(title='% of sequences from taxon<br>aligning with human sequence', row=1,col=1)



    fig.update_layout(legend=dict(
        orientation="v",
        #yanchor="bottom",
        #y=1.02,
        xanchor="left",
        x=1
    ))

    return(fig)

def make_conservation_plot(result):
    from scipy.stats import mannwhitneyu as mu

    # Modern = 1, Ancient = 0
    ancient = result[result['Current Study'] == "Ancient"]["Conservation"].dropna()
    modern = result[result['Current Study'] == "Modern"]["Conservation"].dropna()

    if len(ancient) == 0 or len(modern) == 0:
        mannP = "Undefined"
        mannU = "Undefined"
    else:
        mann = mu(ancient, modern, alternative="greater")
        mannU = f"{mann.statistic:.2f}"
        mannP = f"{mann.pvalue:.2e}"

    hist1 = px.histogram(result, x="Conservation", histnorm='probability', color='Current Study', nbins=100, marginal="rug",barmode='overlay',opacity=0.5, title="Current Study", category_orders={"Current Study": ["Ancient","Modern"]})
    hist1.add_annotation(
        xref="x", yref="y",
        x=0, y=0,
        text=f"Mann-Whitney U test p-value: {mannP}",
        showarrow=False,
        xanchor="left",
        yanchor="top",
    )
    currCountAnc = len(ancient)
    currCountMod = len(modern)
    currMeanAnc = f"{ancient.mean():.2f} ({ancient.sem():.2f})"
    currMeanMod = f"{modern.mean():.2f} ({modern.sem():.2f})"
    currMedAnd = ancient.median()
    currMedMod = modern.median()
    currPvalue = mannP
    currStat = mannU

    # Modern = 1, Ancient = 0
    ancient = result[result['Guo et al (2010)'] == "Ancient"]["Conservation"].dropna()
    modern = result[result['Guo et al (2010)'] == "Modern"]["Conservation"].dropna()

    if len(ancient) == 0 or len(modern) == 0:
        mann_guoP = "Undefined"
        mann_guoU = "Undefined"
    else:
        mann_guo = mu(ancient, modern, alternative="greater")
        mann_guoU = f"{mann_guo.statistic:.2f}"
        mann_guoP = f"{mann_guo.pvalue:.2e}"

    guoCountAnc = len(ancient)
    guoCountMod = len(modern)
    guoMeanAnc = f"{ancient.mean():.2f} ({ancient.sem():.2f})"
    guoMeanMod = f"{modern.mean():.2f} ({modern.sem():.2f})"
    guoMedAnd = ancient.median()
    guoMedMod = modern.median()
    guoPvalue = mann_guoP
    guoStat = mann_guoU


    hist2 = px.histogram(result, x="Conservation", histnorm='probability', color='Guo et al (2010)', nbins=100, marginal="rug",barmode='overlay',opacity=0.5, title="Guo et al (2010)", category_orders={"Guo et al (2010)": ["Ancient","Modern"]})
    hist2.add_annotation(
        xref="x", yref="y",
        x=0, y=0,
        text=f"Mann-Whitney U test p-value: {mann_guoP}",
        showarrow=False,
        xanchor="left",
        yanchor="top",
    )
    
    info = pd.DataFrame([[currCountAnc,currCountMod,currMeanAnc, currMeanMod, currMedAnd, currMedMod, currStat, currPvalue],
                         [guoCountAnc, guoCountMod, guoMeanAnc, guoMeanMod, guoMedAnd, guoMedMod, guoStat, guoPvalue]],
                        index = ["Current Study","Guo et al (2010)"],
                        columns= ["# of sites Ancient", "# sites Modern", "Mean Ancient (SD)","Mean Modern (SD)","Median Ancient","Median Modern","Mann-Whitney U","p-value"])


    return(info, hist1, hist2)


def get_vars(name, entry):
    '''Make table with missense variants from UniProt and hgmd.

    Arguments:
    name -- str, name of the analysis
    entry -- str, uniprot accession number of the reference protein

    '''
    
    # Get variants from UniProt in JSON format
    import requests, sys
    requestURL = f"https://www.ebi.ac.uk/proteins/api/variation/{entry}"
    r = requests.get(requestURL, headers={ "Accept" : "application/json"})
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    responseBody = r.text

    # Parse the JSON response and make a table with the variants
    # Only keep the missense variants
    
    import json
    variants = json.loads(responseBody)
    
    variantsList = []

    for i in variants['features']:
        if i['consequenceType'] == 'missense':
            pos = i['begin']
            
            if 'mutatedType' in i:
                mutation = i['wildType'] + pos + i['mutatedType']
            else:
                mutation = i['wildType'] + pos + 'X'
            consequence = 'missense'
            names = ""
            for j in i['xrefs']:
                if 'name' not in j:
                    names = names + "missing name:"
                else:
                    if 'id' in j:
                        names = names + f"{j['name']}:{j['id']}; "  
            if 'clinicalSignificances' in i:
                if i['clinicalSignificances'][0]['type'] == 'Pathogenic':
                    patho = 'Pathogenic'
                else:
                    patho = 'Benign/Unknown'
            else:
                patho = 'Benign/Unknown'
            text = f"<b>Mutation:</b> {mutation} <br> <b>Type:</b> {consequence} <br> <b>Clinical:</b> {patho}"
            
            # Remove cancer associated variants
            if 'cosmic' in names or 'Cosmic' in names or 'NCI-TCGA' in names:
                next
            else:
                variantsList.append([pos,mutation,consequence,names,patho,text])

    variantsDf = pd.DataFrame(variantsList,columns=['Pos','Mutation','Consequence','Database','Pathogenicity','Text'])

    # Add variants extracted from HGMD
    hgmd = pd.read_excel(f"hgmd/{name}.xlsx")
    # Filter out nonsense variants, "*"; synonymous variants, "="; and variants changing the start codon, "?".
    hgmd = hgmd[~hgmd['Protein Change'].str.contains(r'[=*?]')]

    # Create new columns with the position and the mutation, based on Protein Change column.
    # matching = hgmd["Protein Change"].str.extract(r"p\.(?P<orig>\w)(?P<pos>\d+)(?P<final>\w)")
    # matching.to_csv(f"{name}/matching.csv",index=False)
    # hgmd.to_csv(f"{name}/hgmd.csv",index=False)
    # hgmd["Pos"] = matching['pos']
    # hgmd["Mutation"] = matching['orig'] + "->" + matching['final']
    hgmd["Mutation"] = hgmd["Protein Change"].str.extract(r"p\.(\w\d+\w)")
    hgmd["Pos"] = hgmd["Protein Change"].str.extract(r"p\.\w(\d+)\w")
    hgmd["Consequence"] = "missense"
    hgmd["Database"] = "HGMD:"+hgmd["Accession"]
    hgmd["Pathogenicity"] = "Pathogenic"
    hgmd["Text"] = hgmd['Phenotype']
    
    hgmd = hgmd[["Pos","Mutation","Consequence","Database","Pathogenicity","Text"]]
    
    allVariants = pd.concat([variantsDf,hgmd],ignore_index=True)
    allVariants["Pos"] = allVariants["Pos"].astype(int)
    
    # Remove variants in the start codon, as their effect might be larger than just a missense
    allVariants = allVariants[allVariants['Pos']>1]
    
    # Collapse identical variants
    allVariants = allVariants.drop_duplicates(subset=['Mutation'], keep='last')
    allVariants = allVariants.sort_values(by='Pos')
    
    return(allVariants)

def calculate_scores(name):
    from skbio import TabularMSA, Protein

    # Read the alignment file
    msa = TabularMSA.read(f'{name}/{name}.aln', constructor=Protein)
    
    # Remove all positions that are gaps in the reference sequence, the first sequence in the alignment
    mask = []
    for i in str(msa[0,:]):
        if i == "-":
            mask.append(False)
        else:
            mask.append(True)
    filtered_msa = msa[:, mask]
    
    # Save the filtered alignment
    filtered_msa.write(f'{name}/{name}_filter.aln')
    
    # Calculate and save the conservation scores
    positional_conservation = filtered_msa.conservation(metric='inverse_shannon_uncertainty', degenerate_mode='nan', gap_mode='ignore')
    positional_conservation = pd.DataFrame(positional_conservation,columns=["Conservation"])
    positional_conservation.to_csv(f'{name}/{name}.cons',sep="\t",index=False)
    
    return(positional_conservation)

def align_seqs(name):
    filein = f"{name}/{name}.fasta"
    fileout = f"{name}/{name}.aln"

    import subprocess
    try:
        bstring = f'mafft {filein} > {fileout}' # By default mafft outputs the file in the same order as input
        otuput = subprocess.run(bstring,shell=True,check=True,capture_output=True)
        print(otuput.stdout.decode())
    except:
        print(f"Error with {name}")
    return(0)

def save_fasta(name, table, reference):
    '''Saves a fasta file with sequences of interest, with the reference sequence first, for alignment.
    
    Arguments:
    name -- str, name of the analysis
    table -- pd.DataFrame, table with the proteins from the taxa to analyse that were blast hits after filtering
    reference -- str, uniprot accession number of the reference protein
    
    Returns:
    No return value. The fasta file is saved in the folder with the name of the analysis.
    '''
    
    fileout = open(f"{name}/{name}.fasta","w")
    fastStr = ""
    for i, row in table.iterrows():
        seq = row['Sequence']
        taxon = row['Taxon']
        sname = row['Entry Name']
        org = row['Organism (ID)']

        if sname != reference:
            fastStr = fastStr + f">{taxon}_{org}_{sname}\n{seq}\n"

        else:
            fastStr = f">{taxon}_{org}_{sname}\n{seq}\n" + fastStr

    fileout.write(fastStr)
    fileout.close()
    return(0)

def plot_blast(data,reference,blastTable,orderTax,logfile):
    '''Plot the blast coverage of the reference protein by the proteins in the blast table.
    
    Arguments:
    data -- pd.DataFrame, table with the proteins downloaded from UniProt (filtered)
    reference -- str, uniprot accession number of the reference protein
    blastTable -- pd.DataFrame, table with the filtered blast results
    orderTax -- list, list with the taxonomic groups to be included in the analysis
    logfile -- file, file to log errors and warnings
    '''
    
    # Get length of reference sequence
    reference = data[data["Entry Name"] == reference]['Sequence'].values[0]
    length = len(reference)
    
    # Get the total number of species in each taxonomic group
    temp1 = data.drop_duplicates(subset=['Organism'])
    totals = temp1['Taxon'].value_counts()

    
    # Create a table with the coverage of the reference protein by the proteins in the blast table
    # Initialize dictionaries to store the results
    hugeTable = {} # Coverage of the reference protein by each protein in each taxa in the blast table. 0 absent; 1 present
    results = {} # This table will store the coverage of the reference protein by each taxonomic group
    
    for tax in orderTax:
        hugeTable[tax] = {}
        results[tax] = [0] * (length + 1)
    resultstemp = np.zeros(length + 1)

    current = ""
    totTax = {}
    oldTax = ""
    for row, res in blastTable.iterrows():
        tax = res["taxon"]
        start = res["sstart"]
        end = res["send"]
        species = res["entry"]
        
        # If species initialize the species in the huge table
        if species not in hugeTable[tax]:
            hugeTable[tax][species] = [0] * (length + 1)

        # If this HSP hit is to the same sequence as the previous one, add the coverage to the previous one
        if res[0] == current:
            for i in range(start,end+1):
                resultstemp[i] = 1
        # If this HSP hit is to a different sequence, add the coverage as a row to the results table and start a new coverage row
        else:
            # If this is the first line of the blast table, add the coverage to the current taxonomic group
            if oldTax == "":
                results[tax] += resultstemp # Add a list of 0s and 1s to the results table
            # If not add the coverage to the previous taxonomic group
            else:
                results[oldTax] += resultstemp
            current = res[0]
            oldTax = tax
            totTax[oldTax] = totTax.get(oldTax,0) + 1
            resultstemp = np.zeros(length+1)
            for i in range(start,end+1):
                resultstemp[i] = 1

    res = pd.DataFrame(results)

    taxP = []
    
    for tax in orderTax:
        newTax = tax + " %"
        taxP.append(newTax)
        # Divide the number of hits in each position by the total number of hits in the taxonomic group
        res[newTax] = res[tax] * 100 / totTax[tax]

    # Sanity check
    for tax in orderTax:
        if totTax[tax] != totals[tax]:
            logfile.write(f"\tPlot Error in {tax} {totTax[tax]} {totals[tax]}\n\n")

    # Remove position 0 from the table.
    # res = res.iloc[1: , :]
    # if you want to debug or see more info you add huge table to the return values to get info
    return(pd.Series(totTax).rename("Blast"),res)

def refine_blast(name, reference):
    '''Filter the blast results to keep only the best hit for each organism.
    
    Arguments:
    name -- str, name of the analysis
    reference -- str, uniprot accession number of the reference protein
    
    Returns:
    filterRes = pd.DataFrame, table with the filtered blast results
    '''

    filein = f"{name}/{name}.blast"
    filterRes = []
    
    # Read the results of the blast search
    blastRes = pd.read_csv(filein,sep="\t",header=None)
    blastRes.columns = ["qseqid","sseqid","evalue","length","sstart","send","qstart","qend","qseq","sseq"]
    
    # Add columns with the taxon, organism and entry name, derived from the sequence name in the blast results
    blastRes[["taxon","org","entry"]] = blastRes["qseqid"].str.split("_",n=2,expand=True)
    

    # Find the organism of the reference protein and remove all other proteins from the same organism
    orgref = blastRes['org'][blastRes['entry'] == reference].values[0]
    blastRes = blastRes[~((blastRes['org'] == orgref) & (blastRes['entry'] != reference))]
    
    # Group results by organism and count the number of hits for each organism
    temp = blastRes.groupby("org").count()

    # Loop through the organisms
    for ind, row in temp.iterrows():
        if aa == "GluPro":
            # If only one hit for the organism, keep it
            if row["evalue"] == 1: 
                filterRes.append(blastRes.loc[blastRes["org"] == str(ind)].values.flatten().tolist())
            
            # If more than one hit for the same organism, if it is bacteria or archaea keep all, if it is a vertebrate keep the one with the lowest evalue and longest length
            else:
                temp2 = blastRes.loc[blastRes["org"] == str(ind)].sort_values(by=["evalue","length"],ascending=[True,False])
                if temp2["taxon"].iloc[0] in ["Bacteria","Archaea"]:
                    for i, row in temp2.iterrows():
                        filterRes.append(row.values.flatten().tolist())
                else:
                    # Get the entry name of the best hit
                    best = temp2.iloc[0]["entry"]
                    # Get all the HSPs for the best hit
                    for i, row in temp2[temp2["entry"] == best].iterrows():
                        filterRes.append(row.values.flatten().tolist())
        else:
            # If only one hit for the organism, keep it
            if row["evalue"] == 1: 
                filterRes.append(blastRes.loc[blastRes["org"] == str(ind)].values.flatten().tolist())
            # If more than one hit for the same organism, keep the one with the lowest evalue and longest length
            else:
                # Sort all the HSPs for this organism by evalue and length
                temp2 = blastRes.loc[blastRes["org"] == str(ind)].sort_values(by=["evalue","length"],ascending=[True,False])
                # Get the entry name of the best hit
                best = temp2.iloc[0]["entry"]
                # Get all the HSPs for the best hit
                for i, row in temp2[temp2["entry"] == best].iterrows():
                    filterRes.append(row.values.flatten().tolist())

    return(pd.DataFrame(filterRes, columns = ["qseqid","sseqid","evalue","length","sstart","send","qstart","qend","qseq","sseq","taxon","org","entry"]))

def run_blast(intable, name, reference, logfile):
    '''Run blast search of each protein in the filtered table against the reference protein.
    
    Arguments:
    intable -- pd.DataFrame, table with the filtered proteins
    name -- str, name of the analysis
    reference -- str, uniprot accession number of the reference protein
    
    Returns:
    No return value. Results of the Blast search are saved in the blastfile
    '''

    # Copy the input table to avoid modifying the original table
    table = intable.copy()
    
    # Files used
    qfile = f"{name}/{name}_{reference}.fas"  # name of the file with the query sequence
    dbfile = f"{name}/{name}_others.fas"  # name of the file with the database sequences
    rfile = f"{name}/{name}.blast" # name of the file with the blast results
    
    # Create files with the query and database sequences
    # Extract reference sequence and write to file
    refSeq = table[table["Entry Name"] == reference]['Sequence'].values[0]
    fileout = open(qfile, "w")
    fileout.write(f">ref\n{refSeq}\n")
    fileout.close()
    
    # Extract other sequences and write to file. Sequences names will contain the taxon, organism and entry name.
    count = 0
    fileout = open(dbfile, "w")
    for i, row in table.iterrows():
        count = count + 1
        seq = row['Sequence']
        taxon = row['Taxon']
        sname = row['Entry Name']
        org = row['Organism (ID)']
        if name != reference:
            fileout.write(f">{taxon}_{org}_{sname}\n{seq}\n")
    fileout.close()

    # Run blast
    import subprocess
    try:
        bstring = f'blastp -query {dbfile} -subject {qfile} -outfmt "6 qseqid sseqid evalue length sstart send qstart qend qseq sseq" -out {rfile} -evalue 0.01'
        logfile.write(f"\tRunning {bstring}....")
        output = subprocess.run(bstring,shell=True,check=True,capture_output=True)
        logfile.write("Done!\n\n")
    except:
        print(f"\tError with blast search {name}\n\n")
        return(1)
    return(0)

def filter_and_add_taxa(filein, filtSize, lineages, orderTax):
    '''Filter the table downloaded from UniProt based on the specified criteria.
    
    Arguments:
    filein -- str, name of the file with the table downloaded from UniProt
    filtSize -- int, minimum length of the sequences to be included in the analysis
    lineages -- dict, dictionary with the taxonomic lineage IDs for the taxonomic groups to be included in the analysis
    orderTax -- list, list with the taxonomic groups to be included in the analysis
    
    Returns:
    table -- pd.DataFrame, table with the filtered proteins
    stats -- pd.DataFrame, table with the statistics of the filtering process
    '''
    
    # Read table from file downloaded from UniProt
    table = pd.read_csv(filein,sep="\t")
    
    # Add a column, 'Taxon' with taxonomic information based on the lineage IDs that are defined in the LINEAGES dictionary constant, in the order defined by the ORDERTAX list.
    # Note that when one group is contained in another, the more specific group needs to come after the less specific group in the ORDERTAX list. The 'Taxon' column will be assigned in the order of the ORDERTAX list, so in this case all vertebrates will get a value, which is then rewritten for mammals.
    table['Taxon'] = None
    for tax in orderTax:
        numb = lineages[tax]
        table.loc[table['Taxonomic lineage (Ids)'].str.contains(f" {numb} "), 'Taxon'] = tax
    
    
        
    # Make dataframe with counts of all proteins in each lineage
    initial = table['Taxon'].value_counts().rename("Initial")
    
    # Filter table based on the following criteria:

    # It must have a taxonomic lineage in the Taxon column - remove sequences not belonging to any of the taxonomic groups defined in the LINEAGES dictionary
    table = table[~table['Taxon'].isna()]

    # It must have an EC number in the EC number column - sanity check
    table = table[~table["EC number"].isna()]

    # It must have a gene name in the Gene Names column
    table = table[~table["Gene Names"].isna()]

    # It must have a length greater than the specified filtSize
    table = table[table['Length']>filtSize]
    
    # Remove sequences with X characters
    table = table[~table['Sequence'].str.contains('X')]
    
    # Counts of all proteins in each lineage after filtering
    final = table['Taxon'].value_counts().rename("Final")
    
    # Create dataframe with the filtering statistics
    stats = pd.concat([initial,final],axis=1).fillna(0).astype(int)
    stats['Filtered Out'] = stats['Initial'] - stats['Final']
    stats = stats[['Initial','Filtered Out','Final']]

    return(table,stats)

def downloadUniprot(ec, aa, fileout):
    '''Download all the uniprot entries with the specified EC number.
    
    Arguments:
    ec -- str, the EC number to be searched for in uniprot
    fileout -- str, the name of the file where the results will be saved
    '''
    
    # Download the file from uniprot. The downloaded file will have the following columns:
    # accession, reviewed, id, protein_name, gene_names, organism_name, length, lineage, lineage_ids, organism_id, sequence, ec, annotation_score, go_c, xref_alphafolddb
    import requests
    if aa == "GluPro":
        query = f"https://rest.uniprot.org/uniprotkb/stream?compressed=false&fields=accession%2Creviewed%2Cid%2Cprotein_name%2Cgene_names%2Corganism_name%2Clength%2Clineage%2Clineage_ids%2Corganism_id%2Csequence%2Cec%2Cannotation_score%2Cgo_c%2Cxref_alphafolddb&format=tsv&query=%28ec%3A6.1.1.15%29+OR+%28ec%3A6.1.1.17%29"
    else:
        query = f"https://rest.uniprot.org/uniprotkb/stream?compressed=false&fields=accession%2Creviewed%2Cid%2Cprotein_name%2Cgene_names%2Corganism_name%2Clength%2Clineage%2Clineage_ids%2Corganism_id%2Csequence%2Cec%2Cannotation_score%2Cgo_c%2Cxref_alphafolddb&format=tsv&query=%28%28ec%3A{ec}%29%29"
    
    response = requests.get(query)
    
    open(fileout, 'w').write(response.text)
    return(0)

def find_and_filter_new_domains(seqs,minDomLen = 10):
    i = 1
    domains = []
    domList = ['Ancient'] * len(seqs)
    if 1 in seqs:
        i = seqs.index(1)
        while i < len(seqs):
            try :
                j = seqs.index(0,i)
            except:
                j = len(seqs)
            if j - i > minDomLen:
                domains.append([i+1,j])
                for k in range(i,j):
                    domList[k] = "Modern"
            try:
                i = seqs.index(1,j)
            except:
                i = len(seqs)
    return(domains,domList)


def run_pipeline(name, minSize = 400):
    '''Perform all the analysis for a given enzyme.
    
    Arguments:
    name -- str, name of the analysis, in this case the enzyme name
    minSize -- int, minimum size required for a sequence to be included in the analysis

    Returns: 
    
    '''
        
    # Get the annotations for the given protein
    ec = ECNUMBERS[name][0]
    reference = ECNUMBERS[name][1]
    protein = ECNUMBERS[name][2]

    # Read all annotations for aaRS from the Guo et al (2010) paper
    annotations = pd.read_excel('aaRS_domains.xlsx')
    annotation = []
    for i, row in annotations[annotations['Name'] == name].iterrows():
        annotation.append([row['Start'], row['End']])


    ## All the intermediate files and results from the analysis will be saved in the folder with the name of the analysis. 
    ## If the folder does not exist, it will be created.
    ## If the folder already exists, the analysis will restart using the files already in the folder.
    import os
    if os.path.exists(name):
        print(f"{name} folder already exists. Restarting analysis.")
    else:
        print(f"Creating folder {name}")
        os.mkdir(name)

    # FILE NAMES USED
    # name/name.tsv -- all the proteins with the specified EC number downloaded from uniprot
    # name/name.blast -- blast results
    # name/name.log -- log file with info

    uniprotFile = f"{name}/{name}.tsv"
    uniprotFiltered = f"{name}/{name}_filtered.tsv"
    blastfile = f"{name}/{name}.blast"
    filterblast = f"{name}/{name}_filter.blast"


    logfile = open(f"{name}/{name}.log","a")
    logfile.write(f"\n\n##### {datetime.now()} \nStarting analysis for {name} with EC number {ec} and reference {reference}.\n\n")
     

# Download all the proteins in uniprot with the specified EC number, skip step if the file already exists
    logfile.write("1. Get data from UniProt.\n")
    
    if os.path.exists(uniprotFile):
        logfile.write(f"\t{uniprotFile} already present. Using previously downloaded data.\n\n")
    else:
        logfile.write(f"\tDownloading data from UniProt and saving to {uniprotFile}.....")
        downloadUniprot(ec,aa, uniprotFile)
        logfile.write("Done!\n\n")
    
# Filter the proteins downloaded from UniProt. Uses the LINEAGES and ORDERTAX dictionaries to select only the proteins from the specified taxonomic groups.
    logfile.write("2. Filter data.\n")
    table,stats = filter_and_add_taxa(uniprotFile, minSize, LINEAGES, ORDERTAX)
    table.to_csv(uniprotFiltered,sep="\t",index=False)
    stats.to_string(logfile)
    logfile.write("\n\n")

# Blast all sequences against the reference sequence and save to blast results file, skip step if the file already exists
    logfile.write("3. Running blast search.\n")
    
    if os.path.exists(blastfile):
        logfile.write(f"\tFile {blastfile} already present. Using previous results. Skip this step.\n\n")
    else:
        logfile.write(f"\tRunning blast....")
        run_blast(table, name, reference, logfile)
        logfile.write("Done!\n\n")
    
# Filter blast results
    # Keep a single sequence per organism.
    # If there are multiple sequence hits for the same organism, keep the one with the highest score.
    # Get all HSPs for the best hit in each organism
    logfile.write("4. Filter blast results.\n")
    
    if os.path.exists(filterblast):
        logfile.write("\tFile {filterblast} already present. Using previous results.\n\n")
        result = pd.read_csv(f"{name}/{name}_filter.blast",sep="\t")
    else:
        logfile.write(f"\tFiltering blast and saving filtered table to {filterblast}....")
        result = refine_blast(name,reference)
        result.to_csv(f"{name}/{name}_filter.blast",sep="\t",index=False)
        logfile.write("Done!\n\n")
    
# Create a table with the Blast coverage plot
    logfile.write("5. Making Blast Coverage Plot.\n")
    totTax, resultGraph = plot_blast(table,reference,result,ORDERTAX,logfile)
    stats = pd.concat([stats,totTax],axis=1)
    logfile.write("\tDone!\tStats added to the table.\n")
    stats.to_string(logfile)
    logfile.write("\n\n")

# Add columns to the dataframe with domain annotations from Guo et al 2010 and current study
    # Modern = 1; Ancient = 0
    resultGraph["Current Study Raw"] = np.where((resultGraph["Bacteria %"] < COVERAGE) & (resultGraph["Archaea %"] < COVERAGE), 1, 0)
    
    # Process the newly created columns to remove domains shorter than 10 amino acids and get domain annotation
    domains, domLoc = find_and_filter_new_domains(list(resultGraph["Current Study Raw"]),MINDOMAINLEN)
    resultGraph["Current Study"] = domLoc
    
    tableDomains = pd.DataFrame([str(domains).replace("], [", " | ").replace("[","").replace("]","") if str(domains) != "[]" else "No Domains Found",
                                str(annotation).replace("], [", " | ").replace("[","").replace("]","") if str(annotation) != "[]" else "No Domains Listed"], columns = ["Domains"], index = ["Current Study", "Guo et al (2010)"])
    
    
    # only creating the next column to be able to have the correct title in figure legend for the conservation plot
    resultGraph["Domain Type"] = resultGraph["Current Study"]
    
    # Domains
    protLength = table[table["Entry Name"] == reference]["Length"].values[0]
    New = ["Ancient"] * (protLength + 1) # to account for 0-based indexing
    for annot in annotation:
        if not np.isnan(annot[0]):
            if annot[1] > protLength:
                print(f"Annotation {annot} exceeds protein length {protLength} for {name}.")
                break
            start = int(annot[0])
            end = min(int(annot[1]), protLength)  # Ensure the end does not exceed protLength
            for i in range(start, end+1):
                New[i] = "Modern"
    resultGraph['Guo et al (2010)'] = New

#### CONVERT DOMAIN COLUMNS TO CAETGORICAL WITH TWO LEVELS 'Ancient' and 'Modern' 
    resultGraph['Current Study'] = pd.Categorical(resultGraph['Current Study'], categories=['Ancient','Modern'], ordered=True)
    resultGraph['Guo et al (2010)'] = pd.Categorical(resultGraph['Guo et al (2010)'], categories=['Ancient','Modern'], ordered=True)

# Align the sequences and calculate conservation scores
    logfile.write("6. Align sequences and calculate conservation scores.\n")
    
    # Save input fasta file for alignment
    if os.path.exists(f"{name}/{name}.fasta"):
        logfile.write(f"\tFile: {name}/{name}.fasta already exists. Going to alignment.\n")
    else:
        logfile.write(f"\tSaving fasta {name}/{name}.fasta....")
        # Get all sequences from the sequence table that are in the blast results and in the wanted taxonomic groups
        wantedseqs = table[table['Entry Name'].isin(result['entry'].unique())]
        wantedseqs = wantedseqs[wantedseqs['Taxon'].isin(WANTEDTAX)]
        wantedseqs.to_csv(f"{name}/{name}_Seqs_in_Coverage.tsv",sep="\t",index=False)

        save_fasta(name,wantedseqs,reference)
        logfile.write("Saved fasta file!\n")
    
    # Align the sequences
    if os.path.exists(f"{name}/{name}.aln"):
        logfile.write(f"\tFile:{name}/{name}.aln already exists. No need to align.\n\n")
    else:
        logfile.write("\tAligning....")
        align_seqs(name)
        logfile.write(f"Done! Results saved to {name}/{name}.aln\n\n")
        
    # Calculate conservation scores
    if os.path.exists(f"{name}/{name}.cons"):
        logfile.write(f"\tFile: {name}/{name}.cons already exists. Using previously calculated conservation scores.\n\n")
        conservation = pd.read_csv(f"{name}/{name}.cons",sep="\t")
    else:
        logfile.write(f"\tCalculating conservation scores....")
        conservation = calculate_scores(name)
        logfile.write(f"Done! Results in file {name}/{name}.cons\n\n")
    resultGraph['Conservation'] = [0] + list(conservation['Conservation']) # Add 0 to account for 0 indexing
    resultGraph['Position'] = range(0, len(resultGraph))
    resultGraph = resultGraph.iloc[1:] # Remove first row as the sequence starts at position 1, not 0
    

# Get variants
    # Need to use the accession number not gene name.
    entry = table[table["Entry Name"] == reference]["Entry"].values[0]
    logfile.write("7. Get variants.\n")
    if os.path.exists(f"{name}/{name}.2vars"):
        logfile.write(f"\tFile {name}/{name}.vars already present. Using previously downloaded file.\n\n")
        variants = pd.read_csv(f"{name}/{name}.vars",sep="\t")
    else:
        logfile.write("\tGetting variants....")
        variants = get_vars(name,entry)
        variants.to_csv(f"{name}/{name}.vars",sep="\t",index=False)
        logfile.write("Done! Variants saved to {name}/{name}.vars\n\n")

# Make plots and save results
    logfile.write("8. Making plots and saving results.\n")

    resultGraph.to_csv(f'{name}/{name}_domain_graphs.csv', index=False)
    table.to_csv(f'{name}/{name}_table.csv', index=False)
    stats.to_csv(f'{name}/{name}_stats.csv', index=False)
    
    coverage = make_coverage_plot(resultGraph, annotation,aa)
    coverage.write_image(f"{name}/{name}_coverage.pdf",width=1200,height=400)
    
    # Agreement with Guo et al. (2010)
    contingency_table = pd.crosstab(resultGraph['Current Study'], resultGraph['Guo et al (2010)'], dropna=False)
    contingency = contingency_table.to_numpy()
    agreement_value = (contingency[0,0] + contingency[1,1]) / contingency.sum()

    agreement = px.density_heatmap(
        resultGraph,
        x='Current Study',
        y='Guo et al (2010)',
        text_auto = True,
        color_continuous_scale='blues',
        category_orders={'Current Study': ['Ancient', 'Modern'], 'Guo et al (2010)': ['Ancient', 'Modern']},
        #width=400,
        #height=400
    )
    agreement.add_annotation(
        xref="x", yref="paper",
        x=0, y=1,
        text=f"Agreement: {agreement_value:.3f}",
        showarrow=False,
        xanchor="left",
        yanchor="bottom",
    )
    
    temp = 1
    # Make Current method coverage plot
    uTestsCons, conservation1, conservation2 = make_conservation_plot(resultGraph)
    
# Calculate Stats:
    tables = make_test_table(resultGraph, variants)

    write_report(name, coverage, agreement, uTestsCons, conservation1, conservation2, stats, tables, tableDomains)

    return(stats,table,resultGraph,variants)


aas = ['Ala','Arg','Asn','Asp','Cys','Gln','Gly','His','Ile','Leu','Lys','Met','Ser','Thr','Trp','Tyr','Val', 'PheA', 'PheB', 'GluPro']

for aa in aas:
    run_pipeline(aa, minSize = 0)
    print(f"{aa} done")
