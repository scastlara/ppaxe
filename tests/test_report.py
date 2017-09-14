# -*- coding: utf-8 -*-
'''
Tests for the summary/report of the analyses
'''
from ppaxe import core
from ppaxe import report
from pycorenlp import StanfordCoreNLP
import json

def test_summary_totalcount():
    '''
    Tests totalcount of ProtSummary
    '''
    article_text = """
             MAPK seems to interact with chloroacetate esterase.
             However, MAPK is a better target for peroxydase.
             The thing is, Schmidtea mediterranea is a good model organism because reasons.
             However, cryoglobulin is better.
         """
    article = core.Article(pmid="1234", fulltext=article_text)
    article.extract_sentences()
    for sentence in article.sentences:
        sentence.annotate()
        sentence.get_candidates()
    summary = report.ReportSummary([article])
    summary.protsummary.makesummary()
    assert(summary.protsummary.prot_table['MAPK']['totalcount'] == 2)

def test_summary_intcount():
    '''
    Tests int_count of ProtSummary
    '''
    article_text = """
             MAPK seems to interact with chloroacetate esterase.
             However, MAPK is a better target for peroxydase.
             The thing is, Schmidtea mediterranea is a good model organism because reasons.
             However, cryoglobulin is better.
         """
    article = core.Article(pmid="1234", fulltext=article_text)
    article.extract_sentences()
    for sentence in article.sentences:
        sentence.annotate()
        sentence.get_candidates()
        for candidate in sentence.candidates:
            candidate.predict()
    summary = report.ReportSummary([article])
    summary.protsummary.makesummary()
    assert(summary.protsummary.prot_table['MAPK']['int_count']['left'] == 2)

def test_summary_prottable_tohtml():
    '''
    Tests int_count of ProtSummary
    '''
    article_text = """
             MAPK seems to interact with chloroacetate esterase.
             However, MAPK is a better target for peroxydase.
             The thing is, Schmidtea mediterranea is a good model organism because reasons.
             However, cryoglobulin is better.
         """
    article = core.Article(pmid="1234", fulltext=article_text, year=2015)
    article.extract_sentences()
    for sentence in article.sentences:
        sentence.annotate()
        sentence.get_candidates()
        for candidate in sentence.candidates:
            candidate.predict()
    summary = report.ReportSummary([article])
    summary.protsummary.makesummary()
    thetable = summary.protsummary.table_to_html(sorted_by="int_count")
    reftable = """<table id="prottable">
<thead>
<tr>
<th>Protein</th>
<th>Total count</th>
<th>Int. count</th>
<th>Left count</th>
<th>Right count</th>
</tr>
</thead>
<tbody>
<tr>
<td><a href="http://www.uniprot.org/uniprot/?query=MAPK&sort=score" target="_blank">MAPK</a></td>
<td>2</td>
<td>2</td>
<td>2</td>
<td>0</td>
</tr>
<tr>
<td><a href="http://www.uniprot.org/uniprot/?query=CHLOROACETATE ESTERASE&sort=score" target="_blank">CHLOROACETATE ESTERASE</a></td>
<td>1</td>
<td>1</td>
<td>0</td>
<td>1</td>
</tr>
<tr>
<td><a href="http://www.uniprot.org/uniprot/?query=PEROXYDASE&sort=score" target="_blank">PEROXYDASE</a></td>
<td>1</td>
<td>1</td>
<td>0</td>
<td>1</td>
</tr>
<tr>
<td><a href="http://www.uniprot.org/uniprot/?query=CRYOGLOBULIN&sort=score" target="_blank">CRYOGLOBULIN</a></td>
<td>1</td>
<td>0</td>
<td>0</td>
<td>0</td>
</tr>
</tbody>
</table>"""
    assert(thetable == reftable)

def test_interaction_list():
    '''
    Tests if GraphSummary.makesummary() creates the interaction list correctly
    '''
    article_text = """
             MAPK seems to interact with MAPK4.
             However, Mapk4 interacts directly with MAPK.
             CPP3 is a molecular target of Akt3.
             AKT3 is also known to interact with CPP3.
         """
    article = core.Article(pmid="1234", fulltext=article_text)
    article.extract_sentences()
    for sentence in article.sentences:
        sentence.annotate()
        sentence.get_candidates()
        for candidate in sentence.candidates:
            candidate.predict()
    summary = report.ReportSummary([article])
    summary.graphsummary.makesummary()
    assert(
        len(summary.graphsummary.interactions) == 4 and
        summary.graphsummary.uniqinteractions_count == 2
    )


def test_interaction_table_html():
    '''
    Tests the markdown of the interactions table
    '''
    article_text = """
             MAPK seems to interact with MAPK4.
             However, Mapk4 interacts directly with MAPK.
             CPP3 is a molecular target of Akt3.
             AKT3 is also known to interact with CPP3.
         """
    article = core.Article(pmid="1234", fulltext=article_text, year=2017)
    article.extract_sentences()
    for sentence in article.sentences:
        sentence.annotate()
        sentence.get_candidates()
        for candidate in sentence.candidates:
            candidate.predict()
    summary = report.ReportSummary([article])
    summary.graphsummary.makesummary()
    reftable = """<table id="inttable">
<thead>
<tr>
<th>Confidence</th>
<th>Protein (A)</th>
<th>Protein (B)</th>
<th>Off.symbol (A)</th>
<th>Off.symbol (B)</th>
<th>PMid</th>
<th>Year</th>
<th>Sentence</th>
</tr>
</thead>
<tbody>
<tr>
<td>0.844</td>
<td>Mapk4</td>
<td>MAPK</td>
<td><a href="http://www.uniprot.org/uniprot/?query=MAPK4&sort=score" target="_blank">MAPK4</a></td>
<td><a href="http://www.uniprot.org/uniprot/?query=MAPK&sort=score" target="_blank">MAPK</a></td>
<td><a href="https://www.ncbi.nlm.nih.gov/pubmed/?term=1234" target="_blank">1234</a></td>
<td>2017</td>
<td>However , <span class="prot"> Mapk4 </span> <span class="verb">interacts</span> directly with <span class="prot"> MAPK </span> .</td>
</tr>
<tr>
<td>0.796</td>
<td>CPP3</td>
<td>Akt3</td>
<td><a href="http://www.uniprot.org/uniprot/?query=CPP3&sort=score" target="_blank">CPP3</a></td>
<td><a href="http://www.uniprot.org/uniprot/?query=AKT3&sort=score" target="_blank">AKT3</a></td>
<td><a href="https://www.ncbi.nlm.nih.gov/pubmed/?term=1234" target="_blank">1234</a></td>
<td>2017</td>
<td><span class="prot"> CPP3 </span> <span class="verb">is</span> a molecular target of <span class="prot"> Akt3 </span> .</td>
</tr>
<tr>
<td>0.744</td>
<td>MAPK</td>
<td>MAPK4</td>
<td><a href="http://www.uniprot.org/uniprot/?query=MAPK&sort=score" target="_blank">MAPK</a></td>
<td><a href="http://www.uniprot.org/uniprot/?query=MAPK4&sort=score" target="_blank">MAPK4</a></td>
<td><a href="https://www.ncbi.nlm.nih.gov/pubmed/?term=1234" target="_blank">1234</a></td>
<td>2017</td>
<td><span class="prot"> MAPK </span> <span class="verb">seems</span> to <span class="verb">interact</span> with <span class="prot"> MAPK4 </span> .</td>
</tr>
<tr>
<td>0.714</td>
<td>AKT3</td>
<td>CPP3</td>
<td><a href="http://www.uniprot.org/uniprot/?query=AKT3&sort=score" target="_blank">AKT3</a></td>
<td><a href="http://www.uniprot.org/uniprot/?query=CPP3&sort=score" target="_blank">CPP3</a></td>
<td><a href="https://www.ncbi.nlm.nih.gov/pubmed/?term=1234" target="_blank">1234</a></td>
<td>2017</td>
<td><span class="prot"> AKT3 </span> <span class="verb">is</span> also <span class="verb">known</span> to <span class="verb">interact</span> with <span class="prot"> CPP3 </span> .</td>
</tr>
</tbody>
</table>"""
    htmltable = summary.graphsummary.table_to_html()
    assert(htmltable == reftable)

def test_journal_plots():
    '''
    Tests journal plot
    '''
    article_text1 = """
    Sak binds to p53 , and studies are underway to provide a molecular context for the Sak-p53 interaction.
    By coimmunoprecipitation coupled with mass spectrometry, we demonstrate that AHNAK interacts with dysferlin.
    """
    journal1 = "PLOS ONE"
    article_text2 = """
    Here we show that KLF4 physically interacts with STAT3 upon cytokine-induced phosphorylation of tyrosine 705 ( Y705 ) on STAT3.
    In this study , we report the Grb7 protein interacts with Filamin-a , an actin-crosslinking component of the cell cytoskeleton.
    """
    journal2 = "BMC GENOMICS"

    articles = [
        core.Article(pmid="1234", fulltext=article_text1, journal=journal1, year = 2009),
        core.Article(pmid="4321", fulltext=article_text2, journal=journal2, year = 2016)
    ]

    for article in articles:
        article.predict_interactions()
    summary = report.ReportSummary(articles)
    fig1, fig2, fig3 = summary.journal_plots()
    assert(fig1 and fig2 and fig3)

def test_report_html():
    '''
    Tests journal plot
    '''
    article_text1 = """
    Sak binds to p53 , and studies are underway to provide a molecular context for the Sak-p53 interaction.
    By coimmunoprecipitation coupled with mass spectrometry, we demonstrate that AHNAK interacts with dysferlin.
    """
    journal1 = "PLOS ONE"
    article_text2 = """
    Here we show that KLF4 physically interacts with STAT3 upon cytokine-induced phosphorylation of tyrosine 705 ( Y705 ) on STAT3.
    In this study , we report the Grb7 protein interacts with Filamin-a , an actin-crosslinking component of the cell cytoskeleton.
    """
    journal2 = "BMC GENOMICS"

    articles = [
        core.Article(pmid="1234", fulltext=article_text1, journal=journal1, year = 2009),
        core.Article(pmid="4321", fulltext=article_text2, journal=journal2, year = 2016),
    ]

    for article in articles:
        article.predict_interactions()
    summary = report.ReportSummary(articles)
    summary.make_report("kktest")
