'''
Core classes for ppaxe ppi predictor
'''

import requests
from xml.dom import minidom
from nltk.tag.stanford import StanfordNERTagger
#import json

# CLASSES
# ----------------------------------------------
class Article(object):
    '''
    Article class
        PMID:      PubMed Identifier
        PMCID:     PubMedCentral Identifier
        Abstract:  Full abstract as a string
        Fulltext:  Full text of article as an xml minidom object
        Sentences: List of sentence objects
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
        '''
        if source == "fulltext":
            pass
        else:
            pass

    def count_genes(self):
        '''
        Returns how many times each gene appears.
        Dictionary of gene objects with counts as values
        '''
        pass


# ----------------------------------------------
class Sentence(object):
    '''
    Class for sentences
        Original:  Text of the sentence
        tokenjson: JSON string with the sentence tokenized by SP
        genes:     List of Gene objects in the sentence
    '''
    def __init__(self, original):
        '''
        Original required
        '''
        self.original   = original
        self.tokenjson  = None
        self.genes      = list()

    def tokenize(self):
        '''
        Tokenize (in JSON) the sentence using the Standford Parser
        '''
        pass

    def as_html(self):
        '''
        Writes the tokenized sentence as an HTML
        '''
        pass

# ----------------------------------------------
class Gene(object):
    '''
    Class for genes in sentences
    '''
    def __init__(self, symbol, position):
        self.symbol   = symbol
        self.position = position
        self.synonym  = list()

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
