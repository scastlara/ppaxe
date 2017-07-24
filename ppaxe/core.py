# -*- coding: utf-8 -*-
'''
Core classes for ppaxe ppi predictor
'''

import requests
from xml.dom import minidom
import ner
import re
from pycorenlp import StanfordCoreNLP


NLP = StanfordCoreNLP('http://localhost:9000')

#NER_TAGGER = ner.SocketNER(host='localhost', port=9000)

# CLASSES
# ----------------------------------------------
class PMQuery(object):
    '''
    Class for PubMed queries. Will have Article objects. Will try to
    do only ONE GET request... so that the time to retrieve the articles is reduced
    '''
    def __init__(self, ids, database="PMC"):
        self.ids = ids
        self.database = database
        self.articles = list()
        self.notfound = list()

    def get_articles(self):
        '''
        Retrieves the Fulltext or the abstracts of the specified Articles
        '''
        if self.database == "PMC":
            # Do fulltext query
            params = {
                    'id': ",".join(self.ids),
                    'db': 'pmc',
            }
            req = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi", params=params)
            if req.status_code == 200:
                article_text = minidom.parseString(req.content)
                articles = article_text.getElementsByTagName('article')
                for article in articles:
                    pmid  = article.getElementsByTagName('article-id')[0].firstChild.nodeValue
                    pmcid = article.getElementsByTagName('article-id')[1].firstChild.nodeValue
                    body =  article.getElementsByTagName('body')
                    paragraphs = body[0].getElementsByTagName('p')
                    fulltext = list()
                    for par in paragraphs:
                        fulltext.append(" ".join(t.nodeValue.encode('utf-8') for t in par.childNodes if t.nodeType == t.TEXT_NODE))
                    self.articles.append(Article(pmid=pmid, pmcid=pmcid, fulltext="\n".join(fulltext)))
            else:
                PubMedQueryError("Can't connect to PMC...")

        elif self.database == "PUBMED":
            # Do abstract query
            params = {
                    'id':      ",".join(self.ids),
                    'db':      'pubmed',
            }
            req = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi", params=params)


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
    def __init__(self, pmid, pmcid=None, fulltext=None, abstract=None):
        '''
        Required: PMID
        '''
        self.pmid      = pmid
        self.pmcid     = pmcid
        self.abstract  = abstract
        self.fulltext  = fulltext
        self.sentences = None
        self.genes     = None
        self.annotated = None

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

    def annotate_sentences(self):
        '''
        Uses stanford parser to tag the proteins in the sentence
        '''
        if not self.sentences:
            self.extract_sentences()
        # Call SP
        for sentence in self.sentences:
            self.annotated = NLP.annotate(sentence,properties={'outputFormat':'json'})['sentences'][0]['tokens']
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

class PubMedQueryError(Exception):
    '''
    Exception raised when query to Pubmed or PMC does not return a 200 response.
    '''
    pass

class ConnectionError(Exception):
    '''
    Exception raised when can't connect to online service such as PubMed or PubMedCentral
    '''
    pass
