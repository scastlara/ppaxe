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
import sys
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

def json_to_sentence(json):
    '''
    Takes a sentence in json dict an returns a string
    '''
    pass

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
        for sentence in self.annotated:
            # Get the proteins
            prot_counter = 0
            state        = 0
            prot_list = list()
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
    FEATURES:
        1. Token Distance
        2. Tokens in Sentence
        3. VB token count (only verbs between A and B)
        4. VBD token count (only verbs between A and B)
        5. VBG token count (only verbs between A and B)
        6. VBN token count (only verbs between A and B)
        7. VBP token count (only verbs between A and B)
        8. VBZ token count (only verbs between A and B)
        9. Max verb Score (only verbs between A and B)
        10. Total verb score (sum of scores) (only verbs between A and B)
        11. Distance from A to closest verb (only verbs between A and B)
        12. Distance from A to farthest verb (only verbs between A and B)
        13. Distance from B to closest verb (only verbs between A and B)
        14. Distance from B to farthest verb (only verbs between A and B)
        15. VB token count (all verbs in sentence)
        16. VBD token count (all verbs in sentence)
        17. VBG token count (all verbs in sentence)
        18. VBN token count (all verbs in sentence)
        19. VBP token count (all verbs in sentence)
        20. VBZ token count (all verbs in sentence)
        21. Max verb Score (all verbs in sentence)
        22. Total verb score (sum of scores) (all verbs in sentence)
        23. Distance from A to closest verb (all verbs in sentence)
        24. Distance from A to farthest verb (all verbs in sentence)
        25. Distance from B to closest verb (all verbs in sentence)
        26. Distance from B to farthest verb (all verbs in sentence)
        ...
    '''
    # This will be pre-calculated from the corpora.
    # Now it is like this for testing and developing purposes
    verb_scores = dict({
        "acetylate":1, "acylate":1, "amidate":1, "brominate":1, "biotinylate":1,
        "carboxylate":1, "cysteinylate":1, "farnesylate":1, "formylate":1, "hydroxilate":1,
        "hydroxylate":1, "methylate":1, "demethylate":1, "myristoylate":1, "myristylate":1,
        "palmitoylate":1, "palmitylate":1, "phosphorylate":1, "dephosphorylate":1, "pyruvate":1,
        "nitrosylate":1, "sumoylate":1, "ubiquitinylate":1, "ubiquitinate":1, "dissociate":1,
        "assemble":1, "attach":1, "bind":2, "complex":1, "contact":1, "couple":1, "multimerize":1,
        "multimerise":1, "dimerize":1, "dimerise":1, "interact":4, "precipitate":1, "regulate":1,
        "inhibit":1, "activate":3, "downregulate":2, "down-regulate":2, "suppress":2, "upregulate":2,
        "up-regulate":1, "block":1, "inactivate":1, "induce":1, "modify":1, "overexpress":1, "promote":1,
        "stimulate":1, "substitute":1, "catalyze":1, "cleave":1, "conjugate":1, "disassemble":1,
        "discharge":1, "mediate":1, "modulate":1, "repress":1, "transactivate":1
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
        self.__verb_features("between")
        self.__verb_features("all")
        self.__pos_features("between")
        self.__pos_features("all")
        self.__pos_features("ahead-1-1")
        self.__pos_features("ahead-1-2")
        self.__pos_features("ahead-2-1")
        self.__pos_features("ahead-2-2")
        self.__pos_features("behind-1-1")
        self.__pos_features("behind-1-2")
        self.__pos_features("behind-2-1")
        self.__pos_features("behind-2-2")

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

    def __get_token_pos(self, mode="all"):
        '''
        Returns a string with the token POS annotations for coordinates 'from'-'to'
        '''
        init_coord  = None
        final_coord = None
        if mode == "all":
            '''
            Retrieve POS for all the sentence
            '''
            init_coord  = 0
            final_coord = len(self.prot1.sentence)
        else:
            '''
            Retrieve POS between candidate genes
            '''
            init_coord  = self.between_idxes[0] - 1
            final_coord = self.between_idxes[1] + 1

        pos_str = list()
        subsentence = self.prot1.sentence[init_coord:final_coord]
        for token in subsentence:
            pos_str.append(token['pos'])
        return ",".join(pos_str)

    def __pos_features(self, mode, coord=None):
        '''
        Gets counts of POS tags and adds those features
        '''
        pos_counts = {
            "CC": 0, "LS": 0, "MD": 0,
            "NN": 0, "NNS": 0, "NNP": 0,
            "NNPS": 0, "PDT": 0, "POS": 0,
            "PRP": 0, "PRP$": 0, "RB": 0,
            "RBR": 0, "RBS": 0, "RP": 0,
            "SYM": 0, "TO": 0, "UH": 0,
            "VB": 0, "VBD": 0, "VBG": 0,
            "VBN": 0, "VBP": 0, "VBZ": 0,
            "WDT": 0, "WP": 0, "WP$": 0,
            "WRB": 0, "IN": 0, "DT": 0, ".": 0,
            "JJ" :0, "CD": 0, "": 0, "-LRB-": 0,
            "-RRB-": 0, "JJR": 0, ":": 0, "FW": 0,
            'JJS': 0, "EX": 0, "''": 0
        }

        if mode == "all" or mode == "between":
            for pos in self.__get_token_pos(mode=mode).split(","):
                pos_counts[pos] += 1
        else:
            if mode == "ahead-1-1":
                idx = self.between_idxes[0]
            elif mode == "ahead-1-2":
                idx = self.between_idxes[0] + 1
            elif mode == "behind-1-1":
                idx = self.between_idxes[0] - 2
            elif mode == "behind-1-2":
                idx = self.between_idxes[0] - 3
            elif mode == "ahead-2-1":
                idx = self.between_idxes[1] + 1
            elif mode == "ahead-2-2":
                idx = self.between_idxes[1] + 2
            elif mode == "behind-2-1":
                idx = self.between_idxes[1] - 1
            elif mode == "behind-2-2":
                idx = self.between_idxes[1] - 2

            try:
                pos = self.prot1.sentence[idx]['pos']
                pos_counts[pos] += 1
            except IndexError:
                pass
        for postag in sorted(pos_counts):
            self.features.append(pos_counts[postag])


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


    def __verb_features(self, flag):
        '''
        Computes all the verb features for the candidate.
        It works for all the verbs in the sentence (flag=all)
        and for only the verbs between protein 1 and protein 2 (flag=between)
        '''
        numverbs = dict({
            'VB': 0,  'VBD': 0, 'VBG': 0,
            'VBN': 0, 'VBP': 0, 'VBZ': 0
        })
        maxscore   = 0
        totalscore = 0
        verb_idxes = list()
        someverb_flag = False
        tokens_to_work = list()

        if flag == "between":
            tokens_to_work = self.prot1.sentence[self.between_idxes[0]:self.between_idxes[1]]
        else:
            tokens_to_work = self.prot1.sentence

        for token in tokens_to_work:
            if re.match('VB[DGNPZ]?', token['pos']):
                someverb_flag = True
                numverbs[token['pos']]  += 1
                verb_idxes.append(token['index'])
                if token['lemma'] in InteractionCandidate.verb_scores:
                    totalscore +=  InteractionCandidate.verb_scores[token['lemma']]
                    if InteractionCandidate.verb_scores[token['lemma']] > maxscore:
                        maxscore = InteractionCandidate.verb_scores[token['lemma']]

        # Compute verb distances for the two proteins
        (cl1, far1, cl2, far2) = (0,0,0,0)
        if someverb_flag is True:
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
