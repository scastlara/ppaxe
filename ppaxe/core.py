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
import ppaxe.feature_names as fn
import pkg_resources
import cPickle as pickle
from scipy import sparse
import logging
import markdown

reload(sys)
sys.setdefaultencoding('utf8')

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

def make_md_row(items):
    '''
    Returns a markdown string with the items
    '''
    table_str = ["| ", " | ".join([str(x) for x in items]), " |\n"]
    table_str = "".join(table_str)
    return table_str

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
                    try:
                        pmcid = article.getElementsByTagName('article-id')[1].firstChild.nodeValue
                    except:
                        sys.exit(0)
                        continue
                    body =  article.getElementsByTagName('body')
                    if len(body) == 0:
                        continue
                    self.found.add(pmid)
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
                    'retmode': 'xml'
            }
            req = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi", params=params)
            if req.status_code == 200:
                article_text = minidom.parseString(req.content)
                articles = article_text.getElementsByTagName('PubmedArticle')
                for article in articles:
                    pmid = article.getElementsByTagName('PMID')[0]
                    pmid_text =" ".join(t.nodeValue.encode('utf-8') for t in pmid.childNodes if t.nodeType == t.TEXT_NODE)
                    abstracts = article.getElementsByTagName('AbstractText')
                    abstract_text = list()
                    for abst in abstracts:
                        abstract_text.append(" ".join(t.nodeValue.encode('utf-8') for t in abst.childNodes if t.nodeType == t.TEXT_NODE))
                    abstract_text = "\n".join(abstract_text)
                    if not abstract_text.strip():
                        continue
                    self.found.add(pmid)
                    self.articles.append(Article(pmid=pmid_text, abstract=abstract_text))
                self.notfound = set(self.ids).difference(self.found)
            else:
                PubMedQueryError("Can't connect to PubMed...")
        else:
            logging.error('%s: Incorrect database. Choose "PMC" or "PUBMED"', self.database)

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
        self.sentences  = list()

    def predict_interactions(self, source="fulltext"):
        '''
        Simple wrapper method to avoid calls to multiple methods
        '''
        self.extract_sentences(source=source)
        for sentence in self.sentences:
            sentence.annotate()
            sentence.get_candidates()
            for candidate in sentence.candidates:
                candidate.predict()

    def as_html(self):
        '''
        Writes tokenized sentences as HTML
        '''
        pass

    def extract_sentences(self, mode="split", source="fulltext"):
        '''
        Finds sentence boundaries and saves them as sentence objects
        Does not work very well.
        '''
        text = ""
        if source == "fulltext":
            text = self.fulltext
        else:
            text = self.abstract

        if mode == "no-split":
            # Don't try to separate the sentence.
            # Everything in the text is just one sentence!
            self.sentences.append(Sentence(originaltext=text))
        else:
            caps = "([A-Z])"
            prefixes = "(Mr|Fig|fig|St|Mrs|Ms|Dr)[.]"
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
            for sentence in sentences:
                if not sentence.strip():
                    continue
                self.sentences.append(Sentence(originaltext=sentence))

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



# ----------------------------------------------
class Protein(object):
    '''
    Class for protein in sentences
    '''
    def __init__(self, symbol, positions, sentence):
        #self.disambiguate()
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
        return self.symbol.upper()

    def __str__(self):
        return "%s found in positions %s" % (self.symbol, ":".join([ str(idx) for idx in self.positions ]))


# ----------------------------------------------
class Sentence(object):
    '''
    Class for sentences
    '''
    def __init__(self, originaltext):
        self.originaltext = originaltext
        self.tokens       = list()
        self.tree         = list()
        self.candidates   = list()
        self.proteins     = list()

    def annotate(self):
        '''
        Annotates the genes/proteins in the sentence
        '''
        if not self.originaltext.strip():
            self.tokens = ""
        annotated = json.loads(NLP.annotate(self.originaltext))
        #print(annotated)
        if annotated['sentences']:
            self.tokens = annotated['sentences'][0]['tokens']

    def get_candidates(self):
        '''
        Gets interaction candidates candidates for sentence
        '''
        if not self.tokens:
            self.annotate()
        # Get the proteins
        prot_counter = 0
        state        = 0
        prot_list = list()
        for token in self.tokens:
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
            protein_symbol = " ".join([self.tokens[token_idx - 1]['word'] for token_idx in prot_pos])
            protein = Protein(
                symbol=protein_symbol,
                positions=prot_pos,
                sentence=self
            )
            self.proteins.append(protein)
            prots_in_sentence.append(protein)
        # Create candidates for sentence
        for prot in itertools.combinations(prots_in_sentence, r=2):
            self.candidates.append(InteractionCandidate(prot1=prot[0], prot2=prot[1]))

    def to_html(self):
        '''
        Sentence to HTML string
        '''
        if not self.tokens:
            self.annotate()
        if not self.candidates:
            self.get_candidates()

        state = 0
        html_list = list()
        for token in self.tokens:
            word = token['word']
            word = re.sub("-LRB-", "(", word)
            word = re.sub("-RRB-", ")", word)
            if state == 0:
                if token['ner'] == "O":
                    if re.match('VB[DGNPZ]?', token['pos']):
                        html_list.append('<span class="verb">%s</span>' % word)
                    else:
                        html_list.append(word)
                else:
                    # Init of protein
                    state = 1
                    html_list.append('<span class="prot">')
                    html_list.append(word)
            else:
                if token['ner'] == "O":
                    # End of protein
                    html_list.append("</span>")
                    if re.match('VB[DGNPZ]?', token['pos']):
                        html_list.append('<span class="verb">%s</span>' % word)
                    else:
                        html_list.append(word)
                    state = 0
                else:
                    # Continues protein
                    html_list.append(word)
        return " ".join(html_list)

    def __str__(self):
        return self.originaltext


# ----------------------------------------------
class InteractionCandidate(object):
    '''
    Class for interaction candidates in articles
    FEATURES: in feature_names.py
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

    PRED_FILE = pkg_resources.resource_filename('ppaxe', 'data/RF_scikit.pkl')
    with open(PRED_FILE, 'rb') as f:
        predictor = pickle.load(f)

    def __init__(self, prot1, prot2):
        self.prot1 = prot1
        self.prot2 = prot2
        self.between_idxes = (prot1.positions[-1], prot2.positions[0] - 1)
        self.label    = None
        self.votes    = None
        self.feat_cols = list() # Will store the feature indices
        self.feat_current_col = 0
        self.feat_vals = list() # Will store the feature values for each index
        self.features_sparse = None

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
        self.__prot_count("between")
        self.__prot_count("all")
        self.__keyword_count("between")

        rows = [0 for val in self.feat_vals]
        self.features_sparse = sparse.coo_matrix((self.feat_vals, (rows, self.feat_cols)), shape=(1,178))

    def __token_distance(self):
        '''
        Token distance from protein A to protein B.
        '''
        subsentence = self.prot1.sentence.tokens[self.between_idxes[0]:self.between_idxes[1]]
        if len(subsentence) > 0:
            self.feat_cols.append(self.feat_current_col)
            self.feat_vals.append(len(subsentence))
        self.feat_current_col += 1

    def __total_tokens(self):
        '''
        Total tokens in sentence
        '''
        if len(self.prot1.sentence.tokens) > 0:
            self.feat_cols.append(self.feat_current_col)
            self.feat_vals.append(len(self.prot1.sentence.tokens))
        self.feat_current_col += 1

    def __prot_count(self, mode="all"):
        '''
        Counts the number of times the proteins of the candidate (A and B)
        appear in the sentence (either in the whole sentence [mode="all"] or between
        A and B [mode="between"] )
        '''
        prota_count = 0
        protb_count = 0
        tokens = list()
        if mode == "all":
            tokens = self.prot1.sentence.tokens
        else:
            init_coord  = self.between_idxes[0] - 1
            final_coord = self.between_idxes[1] + 1
            tokens = self.prot1.sentence.tokens[init_coord:final_coord]

        for token in tokens:
            if token['word'] == self.prot1.symbol:
                prota_count += 1
            if token['word'] == self.prot2.symbol:
                protb_count += 1

        # These features are always > 0
        # We always need to add them
        self.feat_cols.append(self.feat_current_col)
        self.feat_current_col += 1
        self.feat_cols.append(self.feat_current_col)
        self.feat_vals.extend([prota_count, protb_count])
        self.feat_current_col += 1

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
            final_coord = len(self.prot1.sentence.tokens)
        else:
            '''
            Retrieve POS between candidate genes
            '''
            init_coord  = self.between_idxes[0] - 1
            final_coord = self.between_idxes[1] + 1

        pos_str = list()
        subsentence = self.prot1.sentence.tokens[init_coord:final_coord]
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
            'JJS': 0, "EX": 0, "''": 0, ",":0
        }

        if mode == "all" or mode == "between":
            for pos in self.__get_token_pos(mode=mode).split(","):
                if pos in pos_counts:
                    pos_counts[pos] += 1
        for postag in sorted(pos_counts):
            if pos_counts[postag] > 0:
                self.feat_cols.append(self.feat_current_col)
                self.feat_vals.append(pos_counts[postag])
            self.feat_current_col += 1

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
            tokens_to_work = self.prot1.sentence.tokens[self.between_idxes[0]:self.between_idxes[1]]
        else:
            tokens_to_work = self.prot1.sentence.tokens

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

        for value in [
            numverbs['VB'],  numverbs['VBD'],
            numverbs['VBG'], numverbs['VBN'],
            numverbs['VBP'], numverbs['VBZ'],
            maxscore, totalscore,
            int(cl1), int(far1),
            int(cl2), int(far2)]:
            if value >= 0:
                self.feat_cols.append(self.feat_current_col)
                self.feat_vals.append(value)
            self.feat_current_col += 1

    def __keyword_count(self, mode="all"):
        '''
        Counts appearance of keywords in sentence
        '''

        keywords = dict({
        'acetylate': 0, 'activate': 0, 'acylate': 0, 'amidate': 0, 'assemble': 0, 'attach': 0,
        'bind': 0, 'biotinylate': 0, 'block': 0, 'brominate': 0, 'carboxylate': 0, 'catalyze': 0,
        'cleave': 0, 'complex': 0, 'conjugate': 0, 'contact': 0, 'couple': 0, 'cysteinylate': 0,
        'demethylate': 0, 'dephosphorylate': 0, 'dimerise': 0, 'dimerize': 0,
        'disassemble': 0, 'discharge': 0, 'dissociate': 0, 'down-regulate': 0, 'downregulate': 0,
        'farnesylate': 0, 'formylate': 0, 'hydroxilate': 0, 'hydroxylate': 0, 'inactivate': 0,
        'induce': 0, 'inhibit': 0, 'interact': 0, 'mediate': 0, 'methylate': 0, 'modify': 0,
        'modulate': 0, 'multimerise': 0, 'multimerize': 0, 'myristoylate': 0, 'myristylate': 0,
        'nitrosylate': 0, 'overexpress': 0, 'palmitoylate': 0, 'palmitylate': 0, 'phosphorylate': 0,
        'precipitate': 0, 'promote': 0, 'pyruvate': 0, 'regulate': 0, 'repress': 0, 'stimulate': 0,
        'substitute': 0, 'sumoylate': 0, 'suppress': 0, 'transactivate': 0, 'ubiquitinate': 0,
        'ubiquitinylate': 0, 'up-regulate': 0, 'upregulate': 0
        })

        tokens = list()
        if mode == "all":
            tokens = self.prot1.sentence.tokens
        else:
            init_coord  = self.between_idxes[0] - 1
            final_coord = self.between_idxes[1] + 1
            tokens = self.prot1.sentence.tokens[init_coord:final_coord]

        for token in tokens:
            if token['lemma'] in keywords:
                keywords[token['lemma']] += 1

        for word, value in sorted(keywords.iteritems()):
            if value > 0:
                self.feat_cols.append(self.feat_current_col)
                self.feat_vals.append(value)
            self.feat_current_col += 1

    def features_todense(self):
        '''
        Returns features as a plain python list. Used for testing.
        '''
        if self.features_sparse is None:
            self.compute_features()
        return self.features_sparse.todense()[0].tolist()[0]

    def predict(self):
        '''
        Computes the votes (prediction) of the candidate
        '''
        if self.features_sparse is None:
            self.compute_features()
        pred = InteractionCandidate.predictor.predict_proba(self.features_sparse)[:,1]
        self.votes = pred
        if pred >= 0.55:
            self.label = True
        else:
            self.label = False

    def to_html(self):
        '''
        Transforms candidate to html with only involved proteins tagged and only
        verbs between proteins tagged
        '''
        init_coord  = self.between_idxes[0]
        final_coord = self.between_idxes[1]
        prot1_coords = [pos - 1 for pos in self.prot1.positions]
        prot2_coords = [pos - 1 for pos in self.prot2.positions]
        html_str = list()
        between = range(init_coord, final_coord)
        for i in range(0, len(self.prot1.sentence.tokens)):
            token = self.prot1.sentence.tokens[i]
            word = token['word']
            word = re.sub("-LRB-", "(", word)
            word = re.sub("-RRB-", ")", word)
            word = re.sub("-LSB-", "[", word)
            word = re.sub("-RSB-", "]", word)
            if i in between:
                if re.match('VB[DGNPZ]?', token['pos']):
                    # Verb in between
                    html_str.append('<span class="verb">%s</span>' % word)
                else:
                    # Normal word in between
                    html_str.append(word)
            elif i in prot1_coords:
                if i == prot1_coords[0]:
                    # Start of protein 1
                    html_str.append('<span class="prot"> %s' % word)
                    if i == prot1_coords[-1]:
                        # Protein of length 1
                        html_str.append('</span>')
                elif i == prot1_coords[-1]:
                    # End of protein of length > 1
                    html_str.append('%s </span>' % word)
                else:
                    # Middle of protein 1
                    html_str.append(word)
            elif i in prot2_coords:
                if i == prot2_coords[0]:
                    # Start of protein 2
                    html_str.append('<span class="prot"> %s' % word)
                    if i == prot2_coords[-1]:
                        # Protein of length 2
                        html_str.append('</span>')
                elif i == prot2_coords[-1]:
                    # End of protein of length > 2
                    html_str.append('%s </span>' % word)
                else:
                    # Middle of protein 2
                    html_str.append(word)
            else:
                # Neither verb nor protein
                html_str.append(word)
        return " ".join(html_str)



    def __str__(self):
        return "[%s] may interact with [%s]" % (self.prot1.symbol, self.prot2.symbol)


class ReportSummary(object):
    '''
    Class for the report summary of the analysis
    '''
    def __init__(self, articles):
        '''
        Initialized either with a PMQuery object or with a list of Article objects
        '''
        try: # Check if articles is a PMQuery
            self.articles = articles.articles
        except AttributeError: # Not a PMQuery
            self.articles = articles
        self.protsummary  = ProteinSummary(self.articles)
        self.graphsummary = GraphSummary(self.articles)
        self.totalarticles = len(self.articles)

    def make_report(self, outfile="report"):
        '''
        Makes all the necessary steps to make the report
        '''
        self.protsummary.makesummary()
        self.graphsummary.makesummary()
        self.write_html(outfile)
        # self.write_markdown(outfile)
        # self.create_pdf(outfile)


    def write_html(self, outfile):
        '''
        Writes a markdown with the report to outfile.
        '''
        outfile = outfile + ".html"
        stylesheet = "https://cdn.rawgit.com/scastlara/ppaxe/51b2e788/ppaxe/data/style.css"
        cytotemplate = "https://cdn.rawgit.com/scastlara/ppaxe/51b2e788/ppaxe/data/cytoscape_template.js"
        with open(outfile, "w") as outf:
            md_str = [
                "# PP-axe Report",
                "## Summary",
                '<div class="sumtable">\n',
                '| | |',
                '|--|--|',
                "| Articles analyzed | %s |" % self.totalarticles,
                "| Proteins found | %s |" % self.protsummary.totalprots,
                "| Interactions retrieved | %s |" % self.graphsummary.numinteractions,
                "| Unique interactions | %s |" % self.graphsummary.uniqinteractions_count,
                '\n</div>',
                "-----",
                "## Interactions",
                '<div class="reptable">\n',
                self.graphsummary.table_to_md(),
                '</div>',
                "-----",
                "## Graph",
                '<div id="cyt"></div>\n'
                "-----",
                "## Proteins",
                '<div class="reptable">\n',
                self.protsummary.table_to_md(),
                '</div>'
            ]
            md_str = "\n".join(md_str)
            extensions = ['extra', 'smarty']
            html_body = markdown.markdown(md_str, extensions=extensions, output_format='html5')
            total_html = [
                '<html>',
                '<meta charset="utf-8" />',
                '<head>',
                '<link rel="stylesheet" type="text/css" href="%s">' % stylesheet,
                '</head>',
                '<body>',
                    '<div id="content">',
                    html_body,
                    '</div>',
                    '<script src="https://code.jquery.com/jquery-2.2.4.min.js"></script>\n',
                    '<script src="%s"></script>\n' % cytotemplate,
                    '''
                    <script>
                        graphelements = %s;
                        cy.load(graphelements);
                        cy.layout( { name: 'cose' } );
                    </script>
                    ''' % self.graphsummary.graph_to_json(),
                '</body>',
                '<html>'
            ]
            outf.write("".join(total_html))

    def create_pdf(self, outfile):
        '''
        Creates a pdf out of a markdown file
        '''
        mdfile = outfile + ".md"
        if not sys.path.isfile(mdfile):
            self.write_html(outfile)
        # Convert Markdown file to pdf

class ProteinSummary(object):
    '''
    Class of the Protein summary for the pdf report.
    Will have:
        - Table with protein ocurrence in sentences, ppi...
        - Mapping to uniprot identifiers when possible.
        - Method to get the most common proteins.
        - Protein p-value ocurrence (compared to all PubMed).
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

    def __md_table_header(self):
        '''
        Returns the header of the markdown protein table as a list
        '''
        colnames = ["Protein","Total count","Int. count", "Left count", "Right count"]
        table_str = make_md_row(colnames)
        table_str = table_str + make_md_row(["-----","-----","-----", "-----", "-----"])
        return table_str

    def table_to_md(self, sorted_by="totalcount", reverse=True):
        '''
        Returns a string with the table in Markdown with proteins sorted by
        sorted_by
        '''
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

        table_str = [self.__md_table_header()]
        for protein in sorted(self.prot_table.items(), reverse=reverse, key=sort_lambda):
            table_str.append(make_md_row([
                protein[0],
                protein[1]['totalcount'],
                str(protein[1]['int_count']['right'] + protein[1]['int_count']['left']) ,
                protein[1]['int_count']['left'],
                protein[1]['int_count']['right']
            ]))
        return "".join(table_str)

    def table_to_html(self, sorted_by="totalcount", reverse=True):
        '''
        Returns an html string with the desired count table by converting
        the markdown table to html
        '''
        mdtbl = self.table_to_md(sorted_by=sorted_by, reverse=reverse)
        extensions = ['extra', 'smarty']
        html = markdown.markdown(mdtbl, extensions=extensions, output_format='html5')
        return html

class GraphSummary(object):
    '''
    Class of the Interactions/Graph summary for the pdf report.
    Will have:
        - The actual interactions.
        - Degree plot.
        - Interactions per Journal name plot.
        - Interactions per year plot.
        - Token Distance plot.
        - Method to write a Cytoscape graph.
    '''
    def __init__(self, articles):
        self.articles = articles
        self.interactions = list()
        self.numinteractions = 0
        self.uniqinteractions = set()
        self.uniqinteractions_count = 0

    def makesummary(self):
        '''
        Makes the summary of the interactions retrieved
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
                                article.pmid
                            ]
                        )
        self.uniqinteractions_count = len(self.uniqinteractions)
        self.interactions     = sorted(self.interactions, key=lambda x: x[0], reverse=True)

    def __md_table_header(self):
        '''
        Returns the header of the markdown interaction table as a list
        '''
        colnames = [
            "Confidence", "Protein (A)","Protein (B)",
            "Off.symbol (A)", "Off.symbol (B)",
            "PMid", "Sentence"
        ]
        table_str = make_md_row(colnames)
        table_str = table_str + make_md_row(["---", "---", "---", "---", "---", "---", "---"])
        return table_str
        #print(table_str)

    def table_to_md(self):
        '''
        Returns a string in markdown with the interactions sorted by votes/confidence
        '''
        table_str = [self.__md_table_header()]
        for interaction in self.interactions:
            table_str.append(make_md_row([
                interaction[0][0],
                interaction[1],
                interaction[3],
                interaction[2],
                interaction[4],
                "[%s](https://www.ncbi.nlm.nih.gov/pubmed/?term=%s)" % (interaction[6], interaction[6]),
                interaction[5]
            ]))
        return "".join(table_str)

    def table_to_html(self):
        '''
        Returns a string in html with the interactions sorted by votes/confidence
        '''
        mdtbl = self.table_to_md()
        extensions = ['extra', 'smarty']
        html = markdown.markdown(mdtbl, extensions=extensions, output_format='html5')
        return html

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
            json_ints.append("{ data: { id: '%s-%s', source: '%s', target: '%s', confidence:'%s', colorEDGE: '#cdbb44', sentence: '%s' }}," % ( interaction[2], interaction[4], interaction[2], interaction[4], interaction[0][0], interaction[5]))
        json_nodes.append("], ")
        json_ints.append("]\n")
        total_json = "{\n" + "\n".join(json_nodes) + "\n".join(json_ints) + "\n}"

        return total_json

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

class ProteinNotFound(Exception):
    '''
    Raised when protein is not in Summary class
    '''
    pass
