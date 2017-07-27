# -*- coding: utf-8 -*-
'''
Core classes for ppaxe ppi predictor
'''

import requests
from xml.dom import minidom
import json
import re
from pycorenlp import StanfordCoreNLP
import itertools
import numpy as np
from bisect import bisect_left, bisect_right
import math

NLP = StanfordCoreNLP('http://localhost:9000')

#NER_TAGGER = ner.SocketNER(host='localhost', port=9000)


# FUNCTIONS
# ----------------------------------------------
def pmid_2_pmc(identifiers):
    '''
    Transforms a list of PubMed Ids to PMC ids
    '''
    pmcids = set()
    params = {
        'ids': ",".join(identifiers),
        'format': 'json'
    }
    req = requests.get("https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/", params=params)
    if req.status_code == 200:
        response = json.loads(req.content)

        for record in response['records']:
            if 'status' in record:
                continue
            pmcids.add(record['pmcid'][3:])
        return list(pmcids)
    else:
        raise PubMedQueryError("Can't convert identifiers through Pubmed idconv tool.")

def take_closest(mylist, mynumber):
    """
    Assumes mylist is sorted. Returns closest value to mynumber.
    If two numbers are equally close, return the smallest number.
    from: https://stackoverflow.com/questions/12141150/from-list-of-integers-get-number-closest-to-a-given-value
    """
    pos = bisect_left(mylist, mynumber)
    if pos == 0:
        return mylist[0]
    if pos == len(mylist):
        return mylist[-1]
    before = mylist[pos - 1]
    after = mylist[pos]
    if after - mynumber < mynumber - before:
        return after
    else:
        return before

def take_farthest(mylist, mynumber):
    """
    Returns the farthest number in a sorted list to mynumber
    """
    distance_first = mylist[0]  - mynumber
    distance_last  = mylist[-1] - mynumber
    if math.fabs(distance_first) > math.fabs(distance_last):
        return mylist[0]
    else:
        return mylist[-1]


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
        self.found    = set()
        self.notfound = set()

    def get_articles(self):
        '''
        Retrieves the Fulltext or the abstracts of the specified Articles
        '''
        if self.database == "PMC":
            # Do fulltext query
            params = {
                    'id': ",".join(pmid_2_pmc(self.ids)),
                    'db': 'pmc',
            }
            req = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi", params=params)
            if req.status_code == 200:
                article_text = minidom.parseString(req.content)
                articles = article_text.getElementsByTagName('article')
                for article in articles:
                    pmid  = article.getElementsByTagName('article-id')[0].firstChild.nodeValue
                    pmcid = article.getElementsByTagName('article-id')[1].firstChild.nodeValue
                    self.found.add(pmid)
                    body =  article.getElementsByTagName('body')
                    if len(body) == 0:
                        continue
                    paragraphs = body[0].getElementsByTagName('p')
                    fulltext = list()
                    for par in paragraphs:
                        fulltext.append(" ".join(t.nodeValue.encode('utf-8') for t in par.childNodes if t.nodeType == t.TEXT_NODE))
                    self.articles.append(Article(pmid=pmid, pmcid=pmcid, fulltext="\n".join(fulltext)))
                self.notfound = set(self.ids).difference(self.found)
            else:
                PubMedQueryError("Can't connect to PMC...")

        elif self.database == "PUBMED":
            # Do abstract query
            params = {
                    'id':      ",".join(self.ids),
                    'db':      'pubmed',
            }
            req = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi", params=params)

    def __iter__(self):
        return iter(self.articles)

    def __getitem__(self, index):
        return self.articles[index]

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
        self.pmid       = pmid
        self.pmcid      = pmcid
        self.abstract   = abstract
        self.fulltext   = fulltext
        self.candidates = list()
        self.sentences  = list()
        self.proteins   = list()
        self.annotated  = list()

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
        fig_letters = "([A-Ka-k])"
        suffixes = "(Inc|Ltd|Jr|Sr|Co)"
        starters = r"(Mr|Mrs|Ms|Dr|He\s|She\s|It\s|They\s|Their\s|Our\s|We\s|But\s|However\s|That\s|This\s|Wherever)"
        acronyms = "([A-Z][.][A-Z][.](?:[A-Z][.])?)"
        websites = "[.](com|net|org|io|gov)"
        species = r"([A-Z])[.] ?([a-z]+)"
        text = " " + text + "  "
        text = text.replace("\n"," ")
        text = re.sub(prefixes,"\\1<prd>",text)
        text = re.sub(websites,"<prd>\\1",text)
        if "Ph.D" in text:
            text = text.replace("Ph.D.","Ph<prd>D<prd>")
        text = re.sub(r"\s" + caps + "[.] "," \\1<prd> ",text)
        text = re.sub(acronyms+" "+starters,"\\1<stop> \\2",text)
        text = re.sub(caps + "[.]" + caps + "[.]" + caps + "[.]","\\1<prd>\\2<prd>\\3<prd>",text)
        text = re.sub(caps + "[.]" + caps + "[.]","\\1<prd>\\2<prd>",text)
        text = re.sub(" "+suffixes+"[.] "+starters," \\1<stop> \\2",text)
        text = re.sub(" "+suffixes+"[.]"," \\1<prd>",text)
        text = re.sub(" " + caps + "[.]"," \\1<prd>",text)
        text = re.sub(digits + caps + "[.]"," \\1<prd>",text)
        text = re.sub(digits + "[.]" + digits,"\\1<prd>\\2",text)
        text = re.sub(digits + "[.]" + fig_letters,"\\1<prd>\\2",text)
        text = re.sub(species, "\\1<prd> \\2", text)
        if "”"    in text:
            text = text.replace(".”","”.")
        if "\""   in text:
            text = text.replace(".\"","\".")
        if "!"    in text:
            text = text.replace("!\"","\"!")
        if "?"    in text:
            text = text.replace("?\"","\"?")
        if "e.g." in text:
            text = text.replace("e.g.","e<prd>g<prd>")
        if "i.e." in text:
            text = text.replace("i.e.","i<prd>e<prd>")
        text = text.replace(".",".<stop>")
        text = text.replace("?","?<stop>")
        text = text.replace("!","!<stop>")
        text = text.replace("<prd>",".")
        sentences = text.split("<stop>")
        #sentences = sentences[:-1]
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
            each annotated sentence:
                - index
                - word
                - lemma
                - pos
                - characterOffsetEnd
                - originalText
                - ner
        '''
        if not self.sentences:
            self.extract_sentences()
        # Call SP
        for sentence in self.sentences:
            if not sentence.strip():
                continue
            annotated = json.loads(NLP.annotate(sentence))
            self.annotated.append(annotated['sentences'][0]['tokens'])
        # Now add annotated genes to self.genes...

    def get_candidates(self):
        '''
        Retrieves proteins and interaction candidates from the annotated Sentences
        '''
        if not self.annotated:
            self.annotate_sentences()
        # Get the proteins
        prot_counter = 0
        state        = 0
        prot_list = list()
        # Get the proteins
        for sentence in self.annotated:
            for token in sentence:
                if token['ner'] == "P":
                    state = 1
                    if len(prot_list) - 1 < prot_counter :
                        prot_list.append(list())
                    prot_list[prot_counter].append(token['index'])
                else:
                    if state == 1:
                        prot_counter += 1
                        state = 0

            # Create protein objects for sentence
            prots_in_sentence = list()
            for prot_pos in prot_list:
                # Must substract 1 to token_idx because Stanford CoreNLP has indexes from 1-n, while
                # python lists go from 0-n
                protein_symbol = " ".join([sentence[token_idx - 1]['word'] for token_idx in prot_pos])
                protein = Protein(
                    symbol=protein_symbol,
                    positions=prot_pos,
                    sentence=sentence
                )
                self.proteins.append(protein)
                prots_in_sentence.append(protein)
            # Create candidates for sentence
            for prot in itertools.combinations(prots_in_sentence, r=2):
                self.candidates.append(InteractionCandidate(prot1=prot[0], prot2=prot[1]))

# ----------------------------------------------
class Protein(object):
    '''
    Class for protein in sentences
    '''
    def __init__(self, symbol, positions, sentence):
        # disambiguate first..?
        self.disambiguate()
        self.symbol = symbol
        self.positions = positions
        self.sentence = sentence
        self.synonym = list()
        self.count = len(positions)

    def disambiguate(self):
        '''
        Method for disambiguating the gene
        (convert it to the approved symbol if possible)
        '''
        pass

    def __str__(self):
        return "%s found in positions %s" % (self.symbol, ":".join([ str(idx) for idx in self.positions ]))


# ----------------------------------------------
class Sentence(object):
    '''
    Class for sentences
    '''

# ----------------------------------------------
class InteractionCandidate(object):
    '''
    Class for interaction candidates in articles
    '''
    # This will be pre-calculated from the corpora.
    # Now it is like this for testing and developing purposes
    verb_scores = dict({
        'interact': 3,
        'activate': 2
    })
    def __init__(self, prot1, prot2):
        self.prot1 = prot1
        self.prot2 = prot2
        self.between_idxes = (prot1.positions[-1], prot2.positions[0] - 1)
        self.features = list() # I will make it a numpy array in the future



    def compute_features(self):
        '''
        Computes all the necessary features to predict if this InteractionCandidate
        is a real interaction
        '''
        # Will call several private methods here
        self.__token_distance()
        self.__total_tokens()
        self.__verb_features()

    def __token_distance(self):
        '''
        Token distance from protein A to protein B.
        '''
        subsentence = self.prot1.sentence[self.between_idxes[0]:self.between_idxes[1]]
        self.features.append(len(subsentence))

    def __total_tokens(self):
        '''
        Total tokens in sentence
        '''
        self.features.append(len(self.prot1.sentence))

    def __verb_distances(self, pidx, vidxes):
        '''
        Computes the three verb distances
        '''
        closest_verb_idx  = take_closest(vidxes, pidx)
        farthest_verb_idx = take_farthest(vidxes, pidx)

        # Now compute token distance from indexes
        closest_distance  = math.fabs(closest_verb_idx -  pidx)
        farthest_distance = math.fabs(farthest_verb_idx -  pidx)
        return (closest_distance, farthest_distance)


    def __verb_features(self):
        '''
        Number of verbs between proteins
        and the two verb scores
        '''
        numverbs = dict({
            'VB': 0,  'VBD': 0, 'VBG': 0,
            'VBN': 0, 'VBP': 0, 'VBZ': 0
        })
        maxscore   = 0
        totalscore = 0
        verb_idxes = list()
        for token in self.prot1.sentence[self.between_idxes[0]:self.between_idxes[1]]:
            if re.match('VB[DGNPZ]?', token['pos']):
                numverbs[token['pos']]  += 1
                verb_idxes.append(token['index'])
                if token['lemma'] in InteractionCandidate.verb_scores:
                    totalscore +=  InteractionCandidate.verb_scores[token['lemma']]
                    if InteractionCandidate.verb_scores[token['lemma']] > maxscore:
                        maxscore = InteractionCandidate.verb_scores[token['lemma']]

        # Compute verb distances for the two proteins
        (cl1, far1) = self.__verb_distances(self.prot1.positions[-1], verb_idxes)
        (cl2, far2) = self.__verb_distances(self.prot2.positions[-1], verb_idxes)

        self.features.extend(
            [
                numverbs['VB'],  numverbs['VBD'],
                numverbs['VBG'], numverbs['VBN'],
                numverbs['VBP'], numverbs['VBZ'],
                maxscore, totalscore,
                int(cl1), int(far1),
                int(cl2), int(far2)
            ]
        )


    def __str__(self):
        return "[%s] may interact with [%s]" % (self.prot1.symbol, self.prot2.symbol)


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
    Exception raised when can't connect to online service such as PubMed
    or PubMedCentral
    '''
    pass
