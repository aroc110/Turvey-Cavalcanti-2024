import pandas as pd
from datetime import datetime
import plotly.graph_objects as go

from plotly.subplots import make_subplots
import numpy as np

def make_html(filename, aas, countSeqs,variants,domains,currentRatios, currentCounts, currentSummary, guoRatios, guoCounts, guoSummary, proportionTestCurrent, proportionTestGuo, wilcoxCur, wilcoxGuo, consWilcoxCurr, consWilcoxGuo, tableDomains):
    fileout = open(filename,"w")
    
    links = "<ul>\n"
    for aa in aas:
        links = links + f'<li><a href="{aa}/{aa}.html">{aa}</a></li>\n'
    links = links + "</ul>\n"
    
    fileout.write(f'''<!DOCTYPE html>
<html>
<head>
    <title>Home</title>
    <style>
        body {{
            font-family: Arial, Helvetica, sans-serif;
        }}
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
        .coverage {{ grid-area: 1 / 1 / 2 / 3; }}
        .confusion {{ grid-area: 2 / 1 / 3 / 2; }}
        .violin {{ grid-area: 2 / 2 / 3 / 3; }}
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
        
        ul {{
            -moz-column-count: 4;
            -moz-column-gap: 20px;
            -webkit-column-count: 4;
            -webkit-column-gap: 20px;
            column-count: 4;
            column-gap: 20px;
            list-style-type: none
        }}
        </style>
</head>
<body>
    <h1>Summary analysis of aminoacyl synthetases.</h1>
    <p>Analysis performed on {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}</p>
    <h2>Links for invidual aaRS's results</h2>
    <div class="mycols">
    {links}
    </div>
    
    <h2>Per Sequence statistics of the run</h2>
    <div class="grid-1">
        <div class="table2">
        <h3>Number of sequences analyzed</h3>
        {countSeqs.to_html(classes="pandas-table")}
        </div>
        <div class="table3">
        <h3>Number of variants</h3>
        {variants.to_html(classes="pandas-table")}
        </div>
    </div>
    
    <h2>Comparison of Modern domains from current study with Guo et al (2010)</h2>
    <div width=80%>
    {tableDomains.to_html(classes="pandas-table",float_format="{:.0f}".format)}
    
    {domains.to_html()}
    <p>Blue bars represent the number of variants in the current study, while red bars represent the number of variants in Guo et al (2010). The agreement between the two studies is shown on the right side of the plot.</p>
    </div>
    
    <h2>Conservation of Ancient and Modern Domains:</h2>
    <p>Wilcoxon signed-rank test testing whether residues in ancient domains are more conserved than those in modern domains.</p>
    <div class="grid-1">
        <div class="table2">
        <h3>Current Study</h3>
        {consWilcoxCurr.to_html(classes="pandas-table",float_format="{:.3e}".format)}
        </div>
        <div class="table3">
        <h3>Guo et al (2010)</h3>
        {consWilcoxGuo.to_html(classes="pandas-table",float_format="{:.3e}".format)}
        </div>

    </div>
    
    
    
    <h2>Enrichment of Pathogenic and Benign/Unknown Variants for Ancient and Modern Domains (pooled data):</h2>
    <p>* Excluding Ala, PheA and PheB proteins</p>
    <div class="grid-1">
        <div class="table2">
        <h3>Current Study</h3>
        {currentSummary.to_html(classes="pandas-table",float_format="{:10.2f}".format)}
        <p>Results of Two Sample Z Test of Proportions testing whether the number of variants per site is larger in modern or ancient domains:</p>
        {proportionTestCurrent.to_html(classes="pandas-table",float_format="{:.3e}".format)}
        </div>
        <div class="table3">
        <h3>Guo et al (2010)</h3>
        {guoSummary.to_html(classes="pandas-table",float_format="{:10.2f}".format)}
        <p>Results of Two Sample Z Test of Proportions testing whether the number of variants per site is larger in modern or ancient domains:</p>
        {proportionTestGuo.to_html(classes="pandas-table",float_format="{:.3e}".format)}
        </div>
    </div>

    <div>
    <h2>Variants per site in Ancient and Modern Domains for Each aaRS:</h2>
    <div class="grid-1">
        <div class="table2">
        <h3>Current Study</h3>
        {currentRatios.to_html(classes="pandas-table",float_format="{:10.2f}".format)}
        <p>Results of Wilcoxon signed-rank test testing whether the number of variants per site is larger in modern or ancient domains (excluding Ala, PheA and PheB):</p>
        {wilcoxCur.to_html(classes="pandas-table",float_format="{:.3e}".format)}
        </div>
        <div class="table3">
        <h3>Guo et al (2010)</h3>
        {guoRatios.to_html(classes="pandas-table",float_format="{:10.2f}".format)}
        <p>Results of Wilcoxon signed-rank test testing whether the number of variants per site is larger in modern or ancient domains (excluding Ala, PheA and PheB):</p>
        {wilcoxGuo.to_html(classes="pandas-table",float_format="{:.3e}".format)}
        </div>
    </div>
    </div>
    
    <h2>Detailed tables:</h2>
    <h3>Current Study</h3>
    {currentCounts.to_html(classes="pandas-table",float_format="{:10.0f}".format)}
    <h3>Guo et al (2010)</h3>
    {guoCounts.to_html(classes="pandas-table",float_format="{:10.0f}".format)}

</body>
</html>''')
    fileout.close()
    return(0)

def make_plot_domains(domsDfs):

    plotsNo = len(domsDfs)

    if plotsNo % 2 == 0:
        lines = plotsNo // 2
    else:
        lines = (plotsNo // 2) + 1

    fig = make_subplots(rows=lines, cols=2, vertical_spacing=0.02, shared_xaxes='columns', specs=[[{"secondary_y": True}]*2]*lines)


    graphNo = 1
    for aa in aas:
        current = domsDfs[aa].copy()
        
        current['Curr'] = pd.Categorical(current['Current Study'], categories=['Ancient','Modern'], ordered=True)
        current['Guo'] = pd.Categorical(current['Guo et al (2010)'], categories=['Ancient','Modern'], ordered=True)
        contingency_table = pd.crosstab(current['Curr'], current['Guo'], dropna=False)
        contingency = contingency_table.to_numpy()
        agreement_value = (contingency[0,0] + contingency[1,1]) * 100 / contingency.sum()

        
        current['Current Study'] = np.where(current['Current Study'] == 'Modern', 1, 0)
        current['Guo et al (2010)'] = np.where(current['Guo et al (2010)'] == 'Modern', 1, 0)

        if graphNo > lines:
            colNo = 2
            rowNo = graphNo - lines
        else:
            colNo = 1
            rowNo = graphNo
        
        size = len(current)
        locy = 'y'+str((rowNo-1)*4+1)
        yax = 'yaxis'+str((graphNo-1)*2+1)
        yax2 = 'yaxis'+str((graphNo-1)*2+1+1)
        #print(yax,yax2)
        locx = 'x'+str(colNo)
        
        fig.add_annotation(
            xref=locx, yref=locy,
            x=-0.1, y=0,
            text=aa,
            showarrow=False,
            xanchor="right",
            yanchor="bottom",
            )
        # Add shapes
        fig.add_shape(type="rect",
                xref=locx, yref=locy,
                x0=1, y0=0, x1=size, y1=2,
                line_color='black', fillcolor='white',
                line_width=2
                
                )


        

        if aa == aas[-1]:    
            fig.add_trace(go.Bar(x=current['Position'], y=current['Current Study'], name='Current Study', marker_color='blue', marker_line_width=0, marker_line_color='blue', width=1, showlegend=True), row=rowNo, col=colNo)
            fig['layout'][yax].update(range=[0,2])
            fig.add_trace(go.Bar(x=current['Position'], y=current['Guo et al (2010)'], name='Guo et al (2010)', marker_color='red', marker_line_width=0, marker_line_color='red', width=1, showlegend=True), row=rowNo, col=colNo,secondary_y=True)
            fig['layout'][yax2].update(range=[-1,1])
        else:
            fig.add_trace(go.Bar(x=current['Position'], y=current['Current Study'], name='Current Study', marker_color='blue', marker_line_width=0, marker_line_color='blue', width=1, showlegend=False), row=rowNo, col=colNo)
            fig['layout'][yax].update(range=[0,2])
            
            fig.add_trace(go.Bar(x=current['Position'], y=current['Guo et al (2010)'], name='Guo et al (2010)', marker_color='red', marker_line_width=0, marker_line_color='red', width=1, showlegend=False), row=rowNo, col=colNo,secondary_y=True)
            fig['layout'][yax2].update(range=[-1,1])
        
        fig.add_annotation(
            xref=locx, yref=locy,
            x=size+10, y=0,
            text=f"{agreement_value:.1f}%",
            showarrow=False,
            xanchor="left",
            yanchor="bottom",
        )

        
        graphNo += 1

    #fig.update_layout(yaxis2=dict(showline=True))
    fig.update_xaxes(showline=False, linewidth=1, linecolor='black', mirror=True)
    fig.update_yaxes(showline=False, linewidth=1, linecolor='black', mirror=True)
    #fig.update_layout(yaxis_range=[0,1])
    fig.update_xaxes(ticks='')
    fig.update_yaxes(ticks='', showticklabels=False)



    fig.update_layout(template="simple_white") #,width=1200,height=400)

    fig.update_layout(legend=dict(
        orientation="h",
        yanchor="top",
        #yref = "paper",
        y=-0.1,
        xanchor="center",
        x=0.5
    ))
    return(fig)

def make_summary_table(table):
    # Sums
    table4 = pd.DataFrame(columns=["All","Ancient","Modern","Ratio Ancient/Modern"])

    tableCounts = table.copy()
    tableCounts.drop(['Sum','Ala','PheA','PheB'],inplace=True)
    tableCounts.loc['Sum'] = tableCounts.sum(axis=0)
    
    # row 1
    r1a = tableCounts.loc['Sum']['All']['Total Residues']
    r1b = tableCounts.loc['Sum']['Ancient']['Total Residues']
    r1c = tableCounts.loc['Sum']['Modern']['Total Residues']
    r1d = r1b/r1c if r1c > 0 else "Undefined"

    # row 2
    r2a = tableCounts.loc['Sum']['All']['Total Variants']
    r2b = tableCounts.loc['Sum']['Ancient']['Total Variants']
    r2c = tableCounts.loc['Sum']['Modern']['Total Variants']
    r2d = r2b/r2c if r2c > 0 else "Undefined"

    # row 3
    r3a = tableCounts.loc['Sum']['All']['Pathogenic Variants']
    r3b = tableCounts.loc['Sum']['Ancient']['Pathogenic Variants']
    r3c = tableCounts.loc['Sum']['Modern']['Pathogenic Variants']
    r3d = r3b/r3c if r3c > 0 else "Undefined"

    # row 4
    r4a = tableCounts.loc['Sum']['All']['Benign/Unknown Variants']
    r4b = tableCounts.loc['Sum']['Ancient']['Benign/Unknown Variants']
    r4c = tableCounts.loc['Sum']['Modern']['Benign/Unknown Variants']
    r4d = r4b/r4c if r4c > 0 else "Undefined"

    table4.loc["Total Residues"] = [r1a,r1b,r1c,r1d]
    table4.loc["Total Variants"] = [r2a,r2b,r2c,r2d]
    table4.loc["Pathogenic Variants"] = [r3a,r3b,r3c,r3d]
    table4.loc["Benign/Unknown Variants"] = [r4a,r4b,r4c,r4d]

    table4['All'] = table4['All'].astype(int)
    table4['Ancient'] = table4['Ancient'].astype(int)
    table4['Modern'] = table4['Modern'].astype(int)

    return(table4)

def make_tables_conservation(aas,allDfs,tableno):
    tableConCurr = pd.DataFrame(columns=["# of sites Ancient", "# sites Modern", "Mean Ancient (SD)","Mean Modern (SD)","Median Ancient","Median Modern","Mann-Whitney U","p-value"])
    tableConGuo = pd.DataFrame(columns=["# of sites Ancient", "# sites Modern", "Mean Ancient (SD)","Mean Modern (SD)","Median Ancient","Median Modern","Mann-Whitney U","p-value"])

    for aa in aas:
        work = allDfs[aa][tableno].copy()
        work.set_index('Unnamed: 0', inplace=True)
        tableConCurr.loc[aa] = work.loc['Current Study']
        tableConGuo.loc[aa] = work.loc['Guo et al (2010)']
    return(tableConCurr,tableConGuo)

def make_counts_table(allDfs,aas):
    countsDf = pd.DataFrame(columns=["Bacteria","Archaea","Vertebrata","Mammals"])
    variantsDf = pd.DataFrame(columns=["Protein Length","Pathogenic Variants","Benign/Unknown Variants"])
    for aa in aas:
        # Protein counts
        work = allDfs[aa][0].copy()
        work.set_index('Unnamed: 0', inplace=True)
        work.loc['Total']= work.sum()
        work= work.T
        countsDf.loc[aa] = work.loc['Blast']
        
        # Variant counts
        work = allDfs[aa][1].copy()
        variantsDf.loc[aa] = work.iloc[0]
        
    return(countsDf,variantsDf)


def make_tables_variants(aas, allDfs, tableno):
    tableStats = pd.DataFrame(columns=["Pathogenic/Site Ancient","Pathogenic/Site Modern","Benign/Unknown/Site Ancient","Benign/Unknown/Site Modern"])
    newAll = pd.DataFrame(columns=["Total Residues","Total Variants","Pathogenic Variants","Benign/Unknown Variants"])
    newAncient = pd.DataFrame(columns=["Total Residues","Total Variants","Pathogenic Variants","Benign/Unknown Variants"])
    newModern = pd.DataFrame(columns=["Total Residues","Total Variants","Pathogenic Variants","Benign/Unknown Variants"])

    for aa in aas:
        work = allDfs[aa][tableno].copy()
        work.set_index('Unnamed: 0', inplace=True)
        work.T
        work
        
        # Ancient vs Modern domains
        # Full Protein
        noResiduesAl = work['All']['Number of Residues']
        variantsAl = work['All']['Number of Variants']
        pathoAl = work['All']['Number of Pathogenic Variants']
        unkAl = work['All']['Number of Benign/Unknown Variants']
        newAll.loc[aa] = [noResiduesAl,variantsAl,pathoAl,unkAl]

        # Ancient
        noResiduesA = work['Ancient Domains']['Number of Residues']
        variantsA = work['Ancient Domains']['Number of Variants']
        pathoA = work['Ancient Domains']['Number of Pathogenic Variants']
        unkA = work['Ancient Domains']['Number of Benign/Unknown Variants']
        newAncient.loc[aa] = [noResiduesA,variantsA,pathoA,unkA]

        # Modern
        noResiduesM = work['Modern Domains']['Number of Residues']
        variantsM = work['Modern Domains']['Number of Variants']
        pathoM = work['Modern Domains']['Number of Pathogenic Variants']
        unkM = work['Modern Domains']['Number of Benign/Unknown Variants']
        newModern.loc[aa] = [noResiduesM,variantsM,pathoM,unkM]


        # stats
        variantsAl = variantsAl/noResiduesAl
        pathoAl = pathoAl/noResiduesAl
        unkAl = unkAl/noResiduesAl
        
        variantsA = variantsA/noResiduesA
        pathoA = pathoA/noResiduesA
        unkA = unkA/noResiduesA
        
        variantsM = variantsM/noResiduesM if noResiduesM > 0 else "No Modern"
        pathoM = pathoM/noResiduesM if noResiduesM > 0 else "No Modern"
        unkM = unkM/noResiduesM if noResiduesM > 0 else "No Modern"
        
        tableStats.loc[aa] = [pathoA,pathoM,unkA,unkM]

    newTb = pd.concat([newAll,newAncient,newModern],axis=1,keys=["All","Ancient","Modern"])
    newTb.loc['Sum'] = newTb.sum(axis=0)

    sumTable = make_summary_table(newTb)
    
    return(tableStats,newTb,sumTable)

def make_table_domains(aas,allDfs,tableno):
    tableDomains = pd.DataFrame(columns=["Current Study","Guo et al (2010)"])

    for aa in aas:
        work = allDfs[aa][tableno].copy()
        work.set_index('Unnamed: 0', inplace=True)
        current = work.loc['Current Study']['Domains']
        guo = work.loc['Guo et al (2010)']['Domains']
        
        tableDomains.loc[aa] = [current.replace(", ","-"),guo.replace(", ","-")]
    return(tableDomains)

aas = ['Ala','Arg','Asn','Asp','Cys','Gln','Gly','His','Ile','Leu','Lys','Met','Ser','Thr','Trp','Tyr','Val','PheA','PheB',"GluPro"]

allDfs = {}
domsDfs = {}
for aa in aas:
    print(aa)
    allDfs[aa] = pd.read_html(f'{aa}/{aa}.html')
    temp = pd.read_csv(f'{aa}/{aa}_domain_graphs.csv')
    domsDfs[aa] = temp[['Current Study','Guo et al (2010)','Position']]

domains = make_plot_domains(domsDfs)
domains.write_image("domains.svg", width=1200, height=400)
domains.write_image("domains.pdf", width=1200, height=400)
countSeqs, variants = make_counts_table(allDfs,aas)
currentRatios, currentCounts,currentSummary = make_tables_variants(aas,allDfs,3)
guoRatios, guoCounts, guoSummary = make_tables_variants(aas,allDfs,4)
consWilcoxCurr, consWilcoxGuo = make_tables_conservation(aas,allDfs,5)
tableDomains = make_table_domains(aas,allDfs,2)

filename = "Summary.html"


### Calculate proportion tests for total number of pathogenic vs ancient variants
from statsmodels.stats.proportion import proportions_ztest
curpatoZ, curpatoPval = proportions_ztest(
    [currentSummary.loc['Pathogenic Variants']['Ancient'],currentSummary.loc['Pathogenic Variants']['Modern']],
    [currentSummary.loc['Total Residues']['Ancient'],currentSummary.loc['Total Residues']['Modern']],
    alternative = 'larger'
)
curunkZ, curunkPval = proportions_ztest(
    [currentSummary.loc['Benign/Unknown Variants']['Ancient'],currentSummary.loc['Benign/Unknown Variants']['Modern']],
    [currentSummary.loc['Total Residues']['Ancient'],currentSummary.loc['Total Residues']['Modern']],
    alternative = 'smaller'
)
statCur = [[curpatoZ,curpatoPval,"More pathogenic variants ancient sites"],[curunkZ,curunkPval,"More Benign/Unknown variants in modern sites"]]
proportionTestCurrent = pd.DataFrame(statCur,columns=['Z-score','p-value',"Hypothesis"],index=['Pathogenic','Benign/Unknown'])



guopatoZ, guopatoPval = proportions_ztest(
    [guoSummary.loc['Pathogenic Variants']['Ancient'],guoSummary.loc['Pathogenic Variants']['Modern']],
    [guoSummary.loc['Total Residues']['Ancient'],guoSummary.loc['Total Residues']['Modern']],
    alternative = 'larger'
)
guounkZ, guounkPval = proportions_ztest(
    [guoSummary.loc['Benign/Unknown Variants']['Ancient'],guoSummary.loc['Benign/Unknown Variants']['Modern']],
    [guoSummary.loc['Total Residues']['Ancient'],guoSummary.loc['Total Residues']['Modern']],
    alternative = 'smaller'
)
statGuo = [[guopatoZ,guopatoPval,"More pathogenic variants ancient sites"],[guounkZ,guounkPval,"More Benign/Unknown variants in modern sites"]]
proportionTestGuo = pd.DataFrame(statGuo,columns=['Z-score','p-value',"Hypothesis"],index=['Pathogenic','Benign/Unknown'])

from scipy.stats import wilcoxon
## Wilcoxon test
curW = currentRatios.copy()
curW = curW.drop(['PheA','PheB','Ala'])
curWilcoxPatho = wilcoxon(curW["Pathogenic/Site Ancient"],curW["Pathogenic/Site Modern"],alternative='greater')
curWilcoxUnk = wilcoxon(curW["Benign/Unknown/Site Ancient"],curW["Benign/Unknown/Site Modern"],alternative='less')

wilcoxCur = pd.DataFrame([[curWilcoxPatho.statistic,curWilcoxPatho.pvalue,"More pathogenic variants ancient sites"],[curWilcoxUnk.statistic,curWilcoxUnk.pvalue,"More Benign/Unknown variants in modern sites"]],columns=['W-statistic','p-value',"Hypothesis"],index=['Pathogenic','Benign/Unknown'])

curG = guoRatios.copy()
curG = curG.drop(['PheA','PheB','Ala'])
guoWilcoxPatho = wilcoxon(curG["Pathogenic/Site Ancient"],curG["Pathogenic/Site Modern"],alternative='greater')
guoWilcoxUnk = wilcoxon(curG["Benign/Unknown/Site Ancient"],curG["Benign/Unknown/Site Modern"],alternative='less')

wilcoxGuo = pd.DataFrame([[guoWilcoxPatho.statistic,guoWilcoxPatho.pvalue,"More pathogenic variants ancient sites"],[guoWilcoxUnk.statistic,guoWilcoxUnk.pvalue,"More Benign/Unknown variants in modern sites"]],columns=['W-statistic','p-value',"Hypothesis"],index=['Pathogenic','Benign/Unknown'])




make_html(filename, aas, countSeqs,variants,domains, currentRatios, currentCounts, currentSummary, guoRatios, guoCounts, guoSummary, proportionTestCurrent, proportionTestGuo, wilcoxCur, wilcoxGuo, consWilcoxCurr, consWilcoxGuo, tableDomains)

