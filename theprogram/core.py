'''
Core classes for theprogram ppi predictor
'''

# ----------------------------------------------
class Article(object):
    '''
    Article class
        PMID:      PubMed Identifier
        PMCID:     PubMedCentral Identifier
        Abstract:  Full abstract as a string
        Fulltext:  Full text of article as a string
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

    def download_fulltext(self, source="PMC"):
        '''
        Adds fulltext to the Article object
        '''
        pass

    def as_html(self):
        '''
        Writes tokenized sentences as HTML
        '''
        pass

    def extract_sentences(self):
        '''
        Finds sentence boundaries and saves them as sentence objects
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
