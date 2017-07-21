'''
Core classes for ppaxe ppi predictor
'''

import requests
from xml.dom import minidom
import ner
import re

NER_TAGGER = ner.SocketNER(host='localhost', port=9000)

# CLASSES
# ----------------------------------------------
class PBQuery(object):
    '''
    Class for PubMed queries. Will have Article objects. Will try to
    do only ONE GET request... so that the time to retrieve the articles is reduced
    '''
    pass

# ----------------------------------------------
class Article(object):
    '''
    Article class
        PMID:      PubMed Identifier
        PMCID:     PubMedCentral Identifier
        Abstract:  Full abstract as a string
        Fulltext:  Full text of article as an xml minidom object
        Sentences: List of sentence strings
    '''
    def __init__(self, pmid, pmcid=None):
        '''
        Required: PMID
        '''
        self.pmid      = pmid
        self.pmcid     = pmcid
        self.abstract  = None
        self.fulltext  = None
        self.sentences = None
        self.genes     = None
        self.tagged    = None

    def download_abstract(self):
        '''
        Adds the abstract
        '''
        params = {
            'id':      self.pmid,
            'db':      'pubmed',
            'retmode': 'text',
            'rettype': 'abstract'
        }
        req = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi", params=params)
        if req.status_code == 200:
            self.abstract = req.content
        else:
            raise TextNotAvailable("Abstract not available for PMID: %s" % self.pmid)

    def download_fulltext(self, source="PMC"):
        '''
        Adds fulltext to the Article object
        '''
        if source == "PMC":
            if not self.pmcid:
                pass
            else:
                params = {
                    'id':      self.pmcid,
                    'db':      'pmc',
                }
                req = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi", params=params)
                if req.status_code == 200:
                    self.fulltext = minidom.parseString(req.content)
                else:
                    raise TextNotAvailable("Abstract not available for PMID: %s" % self.pmid)
        elif source == "epub":
            pass
        else:
            pass

    def as_html(self):
        '''
        Writes tokenized sentences as HTML
        '''
        pass

    def extract_sentences(self, source="fulltext"):
        '''
        Finds sentence boundaries and saves them as sentence objects
        Does not work very well.
        '''
        text = ""
        if source == "fulltext":
            text = self.fulltext
        else:
            text = self.abstract

        if self.tagged:
            text = self.tagged

        caps = "([A-Z])"
        prefixes = "(Mr|St|Mrs|Ms|Dr)[.]"
        digits = "([0-9])"
        suffixes = "(Inc|Ltd|Jr|Sr|Co)"
        starters = "(Mr|Mrs|Ms|Dr|He\s|She\s|It\s|They\s|Their\s|Our\s|We\s|But\s|However\s|That\s|This\s|Wherever)"
        acronyms = "([A-Z][.][A-Z][.](?:[A-Z][.])?)"
        websites = "[.](com|net|org|io|gov)"
        text = " " + text + "  "
        text = text.replace("\n"," ")
        text = re.sub(prefixes,"\\1<prd>",text)
        text = re.sub(websites,"<prd>\\1",text)
        if "Ph.D" in text: text = text.replace("Ph.D.","Ph<prd>D<prd>")
        text = re.sub("\s" + caps + "[.] "," \\1<prd> ",text)
        text = re.sub(acronyms+" "+starters,"\\1<stop> \\2",text)
        text = re.sub(caps + "[.]" + caps + "[.]" + caps + "[.]","\\1<prd>\\2<prd>\\3<prd>",text)
        text = re.sub(caps + "[.]" + caps + "[.]","\\1<prd>\\2<prd>",text)
        text = re.sub(" "+suffixes+"[.] "+starters," \\1<stop> \\2",text)
        text = re.sub(" "+suffixes+"[.]"," \\1<prd>",text)
        text = re.sub(" " + caps + "[.]"," \\1<prd>",text)
        text = re.sub(digits + "[.]" + digits,"\\1<prd>\\2",text)
        if "”"    in text: text = text.replace(".”","”.")
        if "\""   in text: text = text.replace(".\"","\".")
        if "!"    in text: text = text.replace("!\"","\"!")
        if "?"    in text: text = text.replace("?\"","\"?")
        if "e.g." in text: text = text.replace("e.g.","e<prd>g<prd>")
        if "i.e." in text: text = text.replace("i.e.","i<prd>e<prd>")
        text = text.replace(".",".<stop>")
        text = text.replace("?","?<stop>")
        text = text.replace("!","!<stop>")
        text = text.replace("<prd>",".")
        sentences = text.split("<stop>")
        sentences = sentences[:-1]
        sentences = [s.strip() for s in sentences]
        self.sentences = sentences

    def count_genes(self):
        '''
        Returns how many times each gene appears.
        Dictionary of gene objects with counts as values
        '''
        pass

    def make_wordcloud(self):
        '''
        Creates a wordcloud image
        '''
        pass

    def tag_proteins(self, source="fulltext"):
        '''
        Uses stanford parser to tag the proteins in the sentence
        '''
        # Call SP
        if source == "fulltext":
            self.tagged = NER_TAGGER.tag_text(self.fulltext)
        else:
            self.tagged = NER_TAGGER.tag_text(self.abstract)
        # Now add annotated genes to self.genes...


# ----------------------------------------------
class Gene(object):
    '''
    Class for genes in sentences
    '''
    def __init__(self, symbol, positions):
        # disambiguate first..?
        self.disambiguate()
        self.symbol    = symbol
        self.positions = positions
        self.synonym   = list()
        self.count     = len(positions)

    def disambiguate(self):
        '''
        Method for disambiguating the gene (convert it to the approved symbol if possible)
        '''
        pass


# EXCEPTIONS
# ----------------------------------------------
class TextNotAvailable(Exception):
    '''
    Exception raised when the text is not available in PubMed
    '''
    pass

class ConnectionError(Exception):
    '''
    Exception raised when can't connect to online service such as PubMed or PubMedCentral
    '''
    pass


'''
# Get things from minidom article
paragraphs = object.getElementsByTagName('p')
for par in paragraphs:
    print " ".join(t.nodeValue for t in par.childNodes if t.nodeType == t.TEXT_NODE)
    print "\n\n\n++++++\n\n\n"
'''
