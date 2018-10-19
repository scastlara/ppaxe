'''
Classes for report Summary
'''
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import base64
import networkx as nx
from io import BytesIO


# FUNCTIONS
# ----------------------------------------------

def make_html_row(items, header=False):
    '''
    Returns an html row
    '''
    if header is True:
        row_str = ['<tr>', '\n'.join([ "<th>" + str(x) + "</th>" for x in items]), '</tr>']
        return "\n".join(row_str)
    else:
        row_str = ['<tr>', '\n'.join([ "<td>" + str(x) + "</td>" for x in items]), '</tr>']
        return "\n".join(row_str)

def year_extender(years):
    '''
    Fills the missing years in a list of years
    '''
    complete_list = list()
    for i in range(0, len(years)):
        if i == len(years) - 1:
            complete_list.append(int(years[i]))
            break
        diff = int(years[i + 1]) - int(years[i])
        complete_list.extend(range(int(years[i]), int(years[i]) + diff ))
    return [str(year) for year in complete_list]

# CLASSES
# ----------------------------------------------
class ReportSummary(object):
    '''
    Class for the report summary of the analysis.

    Attributes
    ----------
    articles : list or PMQuery, no default
        List of Article objects or PMQuery with Article objects in attribute "articles".

    protsummary : ProteinSummary, no default
        ProteinSummary object of the analysis.

    graphsummary : GraphSummary, no default
        GraphSummary object of the analysis.
    '''
    def __init__(self, articles):
        '''
        Summary of the analysis to create an html or pdf report.

        Parameters
        ----------
        articles : list or PMQuery, required, no default
            List of Article objects or PMQuery with Article objects in attribute "articles".
        '''
        try: # Check if articles is a PMQuery
            self.articles = articles.articles
        except AttributeError: # Not a PMQuery
            self.articles = articles
        self.protsummary  = ProteinSummary(self.articles)
        self.graphsummary = GraphSummary(self.articles)
        self.totalarticles = len(self.articles)
        self.totalsentences = sum([ len(art.sentences) for art in self.articles ])
        self.plots = dict()

    def make_report(self, outfile="report"):
        '''
        Makes all the necessary steps to make the report.

        Parameters
        ----------
        outfile : str, optional, default = "report"
            Filename of the output file. Will append ".html" or ".pdf".
        '''
        self.protsummary.makesummary()
        self.graphsummary.makesummary()
        self.plots['j_int_plot'], self.plots['j_prot_plot'], self.plots['a_year_plot'] = self.journal_plots()
        self.plots['network_plot'] = self.graphsummary.network_plot
        self.write_html(outfile)
        # self.write_markdown(outfile)
        # self.create_pdf(outfile)

    def journal_plots(self):
        '''
        Counts the number of proteins and interactions found in each journal.
        Returns the base65 binary of the journal plots
        '''
        journals = dict()
        years = dict()
        for article in self.articles:
            # Initialize
            if article.journal not in journals:
                journals[article.journal] = dict()
                journals[article.journal]['ints'] = 0
                journals[article.journal]['prots'] = 0
            if article.year not in years:
                years[article.year] = 0
            # Count proteins
            for sentence in article.sentences:
                for proteins in sentence.proteins:
                    journals[article.journal]['prots'] += 1
                # Count interactions
                for candidate in sentence.candidates:
                    if candidate.label is True:
                        journals[article.journal]['ints'] += 1
                        # Count years of interactions in articles
                        years[article.year] += 1
        # Remove journals without ints or prots
        journals_ints  = {k: v for k, v in journals.items() if v['ints'] > 0}
        journals_prots = {k: v for k, v in journals.items() if v['prots'] > 0}
        figure_ints =  self.__make_journal_plots(journals_ints,  mode="ints")
        figure_prots = self.__make_journal_plots(journals_prots, mode="prots")
        figure_years = self.__make_year_plot(years)
        sio_ints  = BytesIO()
        sio_prots = BytesIO()
        sio_years = BytesIO()
        figure_ints.savefig(sio_ints, format="png")
        figure_prots.savefig(sio_prots, format="png")
        figure_years.savefig(sio_years, format="png")
        return sio_ints, sio_prots, sio_years

    def __make_year_plot(self, years):
        '''
        Makes plots of the interactions per year
        '''
        labels = sorted(years.keys())
        labels = year_extender(labels)
        year_n = len(labels)
        ind   = np.arange(year_n)
        width = 0.35 # Arbitrary width of bars
        count  = list()
        for lab in labels:
            if lab in years:
                count.append(years[lab])
            else:
                count.append(0)
        # Make the plot
        fig, axis = plt.subplots()
        axis.bar(ind, count, width, color="#777acd")
        axis.set_xticks(ind)
        axis.set_xticklabels(labels, rotation=75)
        axis.set_title("Interactions retrieved per Year")
        axis.set_ylabel("count")
        fig.tight_layout()
        return fig

    def __make_journal_plots(self, journals, mode):
        '''
        Plot the number of proteins and interactions by journals
        '''
        # Check input
        if mode != "ints" and mode != "prots":
            raise IncorrectPlotName("Can't create plot for %s" % mode)
        # Plot parameters
        labels = sorted(journals.keys(), key=lambda x: journals[x][mode], reverse=False)
        journal_n = len(labels)
        ind   = np.arange(journal_n)
        width = 0.35
        count  = list()
        for lab in labels:
            count.append(journals[lab][mode])
        # Make the plot
        fig, axis = plt.subplots()
        if mode == "ints":
            leglabel = "Interactions"
            color = "#c27f3c"
        else:
            leglabel = "Proteins"
            color = "#50ac72"
        axis.barh(ind, count, width, color=color)
        if journal_n > 20:
            labelsize = 5
        elif journal_n > 10:
            labelsize = 8
        else:
            labelsize = 12
        axis.set_yticks(ind)
        max_label_length = 30 # Need to do this because of matplotlib bug #5456
        labels = [ lab[:max_label_length] + (lab[max_label_length:] and '..') for lab in labels ]
        axis.set_yticklabels(labels,  fontsize=labelsize)
        axis.set_title("%s per Journal" % leglabel)
        axis.set_xlabel("count")
        fig.tight_layout()
        return fig

    def summary_table(self):
        '''
        Creates summary table for html report
        '''
        table_str = ['<table class="summarytable">']
        table_str.append(make_html_row(["Articles Analyzed", self.totalarticles]))
        table_str.append(make_html_row(["Total Sentences", self.totalsentences]))
        table_str.append(make_html_row(["Proteins found", self.protsummary.totalprots]))
        table_str.append(make_html_row(["Interactions retrieved", self.graphsummary.numinteractions]))
        table_str.append(make_html_row(["Unique interactions", self.graphsummary.uniqinteractions_count]))
        table_str.append("</table>")
        return "\n".join(table_str)

    def write_html(self, outfile):
        '''
        Writes a markdown with the report to outfile.

        Parameters
        ----------
        outfile : str, required, no default
            Output filename of the html report. Will append ".html" to it.
        '''
        outfile = outfile + ".html"
        stylesheet     = "https://cdn.jsdelivr.net/gh/scastlara/ppaxe@3035a2c94f724f54bbcb327e480aeb6393989ff0/ppaxe/data/style.css"
        cytotemplate   = "https://cdn.jsdelivr.net/gh/scastlara/ppaxe@51b2e7885ca0553adced1d5d362c6f0315ac12d9/ppaxe/data/cytoscape_template.js"
        datatables_css = "https://cdn.datatables.net/1.10.16/css/jquery.dataTables.min.css"
        datatables_js  = "https://cdn.datatables.net/1.10.16/js/jquery.dataTables.min.js"
        with open(outfile, "w") as outf:
            html_str = [
                '<html>',
                '<head>',
                '<meta charset="UTF-8">',
                '<link rel="stylesheet" type="text/css" href="%s">' % stylesheet,
                '<link rel="stylesheet" type="text/css" href="%s">' % datatables_css,
                '<link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0-beta/css/bootstrap.min.css" integrity="sha384-/Y6pD6FV/Vv2HJnA6t+vslU6fwYXjCFtcEpHbNJ0lyAFsXTsjBbfaDjzALeQsN6M" crossorigin="anonymous">',
                '</head>',
                '<body>',
                    '<div id="content">',
                        '<h1>PP-axe Report</h1>',
                        '<h2>Summary</h1>',
                        self.summary_table(),
                        '<hr>',
                        '<h2>Interactions</h2>',
                        '<div class="reptable">',
                            self.graphsummary.table_to_html(),
                        '</div>',
                        '<hr>',
                        '<h2>Graph</h2>',
                        '<div id="cyt"></div>\n',
                        '<hr>',
                        '<h2>Proteins</h2>',
                        '<div class="reptable">',
                            self.protsummary.table_to_html(),
                        '</div>',
                        '<hr>',
                        '<h2>Plots</h2>',
                        '<div class="plots">',
                        '<img id="j_prot_plot" src="data:image/png;base64,%s"/>' % base64.b64encode(self.plots['j_prot_plot'].getvalue()).strip().decode('utf-8'),
                        '<img id="j_int_plot" src="data:image/png;base64,%s"/>' % base64.b64encode(self.plots['j_int_plot'].getvalue()).strip().decode('utf-8'),
                        '<img id="a_year_plot" src="data:image/png;base64,%s"/>' % base64.b64encode(self.plots['a_year_plot'].getvalue()).strip().decode('utf-8'),
                        '</div>',
                    '</div>',
                    '<script src="https://code.jquery.com/jquery-2.2.4.min.js"></script>\n',
                    '<script src="https://cdnjs.cloudflare.com/ajax/libs/popper.js/1.11.0/umd/popper.min.js" integrity="sha384-b/U6ypiBEHpOf/4+1nzFpr53nxSS+GLCkfwBdFNTxtclqqenISfwAzpKaMNFNmj4" crossorigin="anonymous"></script>',
                    '<script src="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0-beta/js/bootstrap.min.js" integrity="sha384-h0AbiXch4ZDo7tp9hKZ4TsHbi047NrKGLO3SEJAg45jXxnGIfYzk4Si90RDIqNm1" crossorigin="anonymous"></script>',
                    '<script src="%s"></script>\n' % cytotemplate,
                    '<script src="%s"></script>\n' % datatables_js,
                    '''
                    <script>
                        graphelements = %s;
                        cy.load(graphelements);
                        cy.layout( { name: 'cose' } );
                    </script>
                    <script>
                    $(document).ready(function(){
                        $('#inttable').DataTable({
                            "order": [[ 0, "desc" ]]
                        });
                        $('#prottable').DataTable({
                            "order": [[ 1, "desc" ]]
                        });
                    });
                    </script>
                    ''' % self.graphsummary.graph_to_json(),
                '</body>',
                '<html>'
            ]
            outf.write("\n".join(html_str))

    def create_pdf(self, outfile):
        '''
        Creates a pdf out of a markdown file.

        Parameters
        ----------
        outfile : str, required, no default
            Output filename of the pdf report. Will append ".pdf" to it.
        '''
        mdfile = outfile + ".md"
        if not sys.path.isfile(mdfile):
            self.write_html(outfile)
        # Convert Markdown file to pdf

class ProteinSummary(object):
    '''
    Class of the Protein summary for the pdf report.

    Attributes
    ----------
    articles : list, no default
        List of Article objects with Article objects in attribute "articles".

    prot_table : dict, no default.
        Dictionary of dictionary with information about protein counts in articles.
        keys:
            symbol: symbol of the protein
                'totalcount' : total number of ocurrencies of protein.
                'int_count'
                    'left'  : Ocurrencies of protein on left hand side of interaction.
                    'right' : Ocurrencies of protein on right hand side of interaction.
    '''
    def __init__(self, articles):
        self.articles = articles
        self.prot_table = dict()
        self.totalprots = 0

    def makesummary(self):
        '''
        Makes the summary of the proteins found using the NER
        '''
        for article in self.articles:
            for sentence in article.sentences:
                for prot in sentence.proteins:
                    symbol = prot.disambiguate()
                    if symbol not in self.prot_table:
                        self.totalprots += 1
                        self.prot_table[symbol] = dict()
                        self.prot_table[symbol]['totalcount'] = 0
                        self.prot_table[symbol]['art_count']  = dict()
                        self.prot_table[symbol]['int_count']  = dict()
                        self.prot_table[symbol]['int_count']['left']  = 0
                        self.prot_table[symbol]['int_count']['right'] = 0
                    self.prot_table[symbol]['totalcount'] += 1
                    if article.pmid not in self.prot_table[symbol]['art_count']:
                        self.prot_table[symbol]['art_count'][article.pmid] = 0
                    self.prot_table[symbol]['art_count'][article.pmid] += 1
                for candidate in sentence.candidates:
                    prot1 = candidate.prot1.disambiguate()
                    prot2 = candidate.prot2.disambiguate()
                    if candidate.label is True:
                        self.prot_table[prot1]['int_count']['left'] += 1
                        self.prot_table[prot2]['int_count']['right'] += 1

    def table_to_html(self, sorted_by="totalcount", reverse=True):
        '''
        Returns an html string with the desired protein/gene count table.

        Parameters
        ----------
        sorted_by : str, optional, default = "totalcount"
            Sort table by total number of ocurrences of protein in sentences (sorted_by="totalcount"),
            by total number of ocurrences in interactions (sorted_by="int_count"), by ocurrences in
            left hand side of interaction (sorted_by="left") or righ hand side (sorted_by="right").

        reverse : bool, optional, default = True
            Sort proteins in reverse order (from bigger to smaller) according to the sorted_by rule if True.
            Reverse (smaller to bigger) if False.
        '''
        # HEADER
        colnames = ["Protein","Total count","Int. count", "Left count", "Right count"]
        table_str = ['<table id="prottable">']
        table_str.append("<thead>")
        table_str.append(make_html_row(colnames, header=True))
        table_str.append("</thead>")

        # Sorted by
        if sorted_by == "totalcount":
            sort_lambda = lambda x: (x[1]['totalcount'], x[1]['int_count']['right'] + x[1]['int_count']['left'])
        elif sorted_by == "int_count":
            sort_lambda = lambda x: x[1]['int_count']['right'] + x[1]['int_count']['left']
        elif sorted_by == "left":
            sort_lambda = lambda x: x[1]['int_count']['right']
        elif sorted_by == "right":
            sort_lambda = lambda x: x[1]['int_count']['right']
        else:
            raise KeyError("Can't sort by %s. Only 'totalcount', 'int_count', 'left' or 'right'.")

        table_str.append("<tbody>")
        for protein in sorted(self.prot_table.items(), reverse=reverse, key=sort_lambda):
            table_str.append(make_html_row([
                '<a href="http://www.uniprot.org/uniprot/?query=%s&sort=score" target="_blank">%s</a>' % (protein[0], protein[0]),
                protein[1]['totalcount'],
                str(protein[1]['int_count']['right'] + protein[1]['int_count']['left']) ,
                protein[1]['int_count']['left'],
                protein[1]['int_count']['right']
            ]))
        table_str.append("</tbody>")
        table_str.append("</table>")
        html = "\n".join(table_str)
        return html

class GraphSummary(object):
    '''
    Class of the Interactions/Graph summary for the pdf report.

    Attributes
    ----------
    articles : list, no default
        List of Article objects with Article objects in attribute "articles".

    interactions : list, no default
        List of lists with interactions in articles.
        elements:
            [
                [
                    votes,
                    prot1.symbol,
                    prot1.disambiguate(),
                    prot2.symbol,
                    prot2.disambiguate(),
                    candidate.to_html(),
                    article.pmid,
                    article.year
                ],
                ...
            ]

    numinteractions : int, no default
        Number of interactions in articles.

    uniqinteractions : set, no default
        Set with symbols of interactions in articles to remove redundant interactions.

    uniqinteractions_count : int, no default
        Number of unique interactions in articles.
    '''
    def __init__(self, articles):
        '''
        Parameters
        ----------
        articles : list, required, no default
            List of Article objects.
        '''
        self.articles = articles
        self.interactions = list()
        self.numinteractions = 0
        self.uniqinteractions = set()
        self.uniqinteractions_count = 0
        self.network_plot = None

    def makesummary(self):
        '''
        Makes the summary of the interactions retrieved.
        '''
        for article in self.articles:
            for sentence in article.sentences:
                for candidate in sentence.candidates:
                    if candidate.label is True:
                        self.numinteractions += 1
                        self.uniqinteractions.add(
                            tuple(sorted([candidate.prot1.disambiguate(), candidate.prot2.disambiguate()]))
                        )
                        self.interactions.append(
                            [
                                candidate.votes,
                                candidate.prot1.symbol,
                                candidate.prot1.disambiguate(),
                                candidate.prot2.symbol,
                                candidate.prot2.disambiguate(),
                                candidate.to_html(),
                                article.pmid,
                                article.year
                            ]
                        )
        self.uniqinteractions_count = len(self.uniqinteractions)
        self.interactions     = sorted(self.interactions, key=lambda x: x[0], reverse=True)


    def make_networkx_plot(self):
        '''
        Saves network plot for PDF output
        '''
        G = nx.Graph()
        ints = list()
        for interaction in self.interactions:
            ints.append((interaction[2], interaction[4]))
        G.add_edges_from(ints)
        options = {
            'node_color': '#FFB600',
            'node_size': 500,
            'with_labels': True,
            'width': 1,

        }
        fig, axis = plt.subplots()
        nx.draw_spring(G, **options)
        sio_net = BytesIO()
        fig.savefig(sio_net, format="png")
        return sio_net


    def table_to_html(self):
        '''
        Returns a string in html with the interactions sorted by votes/confidence
        '''
        colnames = [
            "Confidence", "Protein (A)","Protein (B)",
            "Off.symbol (A)", "Off.symbol (B)",
            "PMid", "Year", "Sentence"
        ]
        table_str = ['<table id="inttable">']
        table_str.append("<thead>")
        table_str.append(make_html_row(colnames, header=True))
        table_str.append("</thead>")

        table_str.append("<tbody>")
        for interaction in self.interactions:
            table_str.append(make_html_row([
                interaction[0],
                interaction[1],
                interaction[3],
                '<a href="http://www.uniprot.org/uniprot/?query=%s&sort=score" target="_blank">%s</a>' % (interaction[2], interaction[2]),
                '<a href="http://www.uniprot.org/uniprot/?query=%s&sort=score" target="_blank">%s</a>' % (interaction[4], interaction[4]),
                '<a href="https://www.ncbi.nlm.nih.gov/pubmed/?term=%s" target="_blank">%s</a>' % (interaction[6], interaction[6]),
                interaction[7],
                interaction[5]
            ]))
        table_str.append("</tbody>")
        table_str.append("</table>")
        return "\n".join(table_str)

    def graph_to_json(self):
        '''
        Returns a json string with the graph prepared for cytoscape
        '''
        json_nodes = list(["nodes: ["])
        json_ints  = list(["edges: ["])
        total_json = list()
        for interaction in self.interactions:
            json_nodes.append("{ data: { id: '%s', name: '%s', colorNODE: '#4b849d' } }," % (interaction[2], interaction[2]))
            json_nodes.append("{ data: { id: '%s', name: '%s', colorNODE: '#4b849d' } }," % (interaction[4], interaction[4]))
            json_ints.append("{ data: { id: '%s-%s', source: '%s', target: '%s', confidence:'%s', colorEDGE: '#cdbb44' }}," % ( interaction[2], interaction[4], interaction[2], interaction[4], interaction[0]))
        json_nodes.append("], ")
        json_ints.append("]\n")
        total_json = "{\n" + "\n".join(json_nodes) + "\n".join(json_ints) + "\n}"

        return total_json


class IncorrectPlotName(Exception):
    '''
    Raised when attempting to create a plot that does not exist
    '''
    pass
