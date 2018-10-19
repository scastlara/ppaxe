# -*- coding: utf-8 -*-
'''
Core classes for ppaxe ppi predictor
'''

import requests
from xml.dom import minidom
import json
import re
import time
from pycorenlp import StanfordCoreNLP
import itertools
from bisect import bisect_left
import math
import sys
import pkg_resources
from scipy import sparse
import logging
import warnings
warnings.filterwarnings("ignore", category=UserWarning)

try:
    # For python 2.7
    import cPickle as pickle
    from HTMLParser import HTMLParser
    reload(sys)
    sys.setdefaultencoding('utf8')
except:
    # For python 3
    import html
    import _pickle as pickle
    from importlib import reload


NLP = StanfordCoreNLP('http://localhost:9000')

# FUNCTIONS
# ----------------------------------------------
def pmid_2_pmc(identifiers):
    '''
    Transforms a list of PubMed Ids to PMC ids
    '''
    pmcids = set()
    maxidents = 200

    for subset in [identifiers[x:x+maxidents] for x in range(0, len(identifiers),maxidents)]:
        params = {
            'ids': ",".join(subset),
            'format': 'json'
        }
        req = requests.get("https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/", params=params)
        if req.status_code == 200:
            response = json.loads(req.content.decode('latin1'))
            for record in response['records']:
                if 'status' in record:
                    continue
                pmcids.add(record['pmcid'][3:])
        else:
            raise PubMedQueryError("Can't convert identifiers through Pubmed idconv tool.")
        time.sleep(3)
    return list(pmcids)

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

def minidom_to_text(minidom):
    '''
    Takes a minidom object and returns the text contained in it without the tags
    '''
    return " ".join(t.nodeValue for t in minidom.childNodes if t.nodeType == t.TEXT_NODE)

# CLASSES
# ----------------------------------------------
class PMQuery(object):
    '''
    Class for PubMed queries. Will have Article objects. Will try to
    do only ONE GET request... so that the time to retrieve the articles is reduced.

    Attributes
    ----------
    ids : list, no default
        List of PubMed identifiers to query.

    database : str, default = "PMC"
        Database to download the articles or abstracts from. PMC or PUBMED.

    articles : list, no default
        List of downloaded Article objects.

    found : set, no default
        PubMed identifiers of the articles found in database.

    notfound : set, no default
        PubMed identifiers of the articles not found in database.

    '''
    def __init__(self, ids, database="PMC"):
        '''
        Parameters
        ----------
        ids : list, required, no default
            List of PubMed identifiers. Required. No default.

        database: str, optional, default = "PMC"
            Database to download the articles or the Abstracts. Can be PMC or PUBMED.
        '''
        self.ids = ids
        self.database = database
        self.articles = list()
        self.found    = set()
        self.notfound = set()

    def __get_pmc(self, req):
        '''
        Gets PMC article request and fills articles attribute

        Parameters
        ----------
        req : requests.models.Response, required, no default
            response object to pubmedCentral
        '''
        if req.status_code == 200:
            article_text = minidom.parseString(req.content)
            articles = article_text.getElementsByTagName('article')
            for article in articles:
                pmid    = article.getElementsByTagName('article-id')[0].firstChild.nodeValue
                journal = article.getElementsByTagName('journal-id')[0].firstChild.nodeValue
                year = article.getElementsByTagName('year')[0].firstChild.nodeValue
                try:
                    pmcid = article.getElementsByTagName('article-id')[1].firstChild.nodeValue
                except:
                    continue
                body =  article.getElementsByTagName('body')
                if len(body) == 0:
                    continue
                self.found.add(pmid)
                paragraphs = body[0].getElementsByTagName('p')
                fulltext = list()
                for par in paragraphs:
                    fulltext.append(minidom_to_text(par))
                self.articles.append(Article(pmid=pmid, pmcid=pmcid, journal=journal, year=year, fulltext="\n".join(fulltext)))
            self.notfound = set(self.ids).difference(self.found)
        else:
            PubMedQueryError("Can't connect to PMC...")

    def __get_pubmed(self, req):
        '''
        Gets PUBMED article request and fills article attribute.

        Parameters
        ----------
        req : requests.models.Response, required, no default
            response object to pubmed
        '''
        if req.status_code == 200:
            article_text = minidom.parseString(req.content)
            articles = article_text.getElementsByTagName('PubmedArticle')
            for article in articles:
                pmid = article.getElementsByTagName('PMID')[0]
                pmid_text = minidom_to_text(pmid)
                journal = article.getElementsByTagName('Journal')[0].getElementsByTagName('Title')[0]
                journal_text = minidom_to_text(journal)
                year = article.getElementsByTagName('Year')[0].firstChild.nodeValue
                abstracts = article.getElementsByTagName('AbstractText')
                abstract_text = list()
                for abst in abstracts:
                    abstract_text.append(minidom_to_text(abst))
                abstract_text = "\n".join(abstract_text)
                if not abstract_text.strip():
                    continue
                self.found.add(pmid_text)
                self.articles.append(Article(pmid=pmid_text, journal=journal_text, year=year, abstract=abstract_text))
            self.notfound = set(self.ids).difference(self.found)
        else:
            PubMedQueryError("Can't connect to PubMed...")

    def get_articles(self):
        '''
        Retrieves the Fulltext or the abstracts of the specified Articles
        '''
        maxidents = 200 # max number of articles per GET request

        for subset in [self.ids[x:x+maxidents] for x in range(0, len(self.ids), maxidents)]:
            if self.database == "PMC":
                # Do fulltext query

                params = {
                    'id': ",".join(pmid_2_pmc(subset)),
                    'db': 'pmc',
                }
                req = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi", params=params)
                self.__get_pmc(req)
            elif self.database == "PUBMED":
                # Do abstract query
                params = {
                    'id':      ",".join(subset),
                    'db':      'pubmed',
                    'retmode': 'xml'
                }
                req = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi", params=params)
                self.__get_pubmed(req)
            else:
                logging.error('%s: Incorrect database. Choose "PMC" or "PUBMED"', self.database)

    def __iter__(self):
        return iter(self.articles)

    def __getitem__(self, index):
        return self.articles[index]

# ----------------------------------------------
class Article(object):
    '''
    Article class.

    Attributes
    ----------
    pmid : str, no default
        PubMed identifier of the article.

    pmcid : str, no default
        PubMedCentral identifier of the article.

    fulltext : str, no default
        Whole text of the article.

    abstract : str, no default
        Abstract of the article.

    journal : str, no default
        Journal id of the article.

    year : str, no default
        Year of publication

    sentences : list, no default
        List of Sentence objects in article (fulltext or abstract).
    '''
    def __init__(self, pmid, pmcid=None, journal=None, year=None, fulltext=None, abstract=None):
        '''
        Parameters
        ----------
        pmid : str, required, no default
            PubMed identifier of the article.

        pmcid : str, optional, no default
            PubMedCentral identifier of the article.

        fulltext : str, optional, no default
            Whole text of the article.

        abstract : str, optional, no default
            Abstract of the article.
        '''
        self.pmid       = pmid
        self.pmcid      = pmcid
        self.journal    = journal
        self.year       = year
        self.abstract   = abstract
        self.fulltext   = fulltext
        self.sentences  = list()

    def predict_interactions(self, source="fulltext", only_dict=False):
        '''
        Simple wrapper method to avoid calls to multiple methods.

        Parameters
        ----------
        source : str, optional, default = "fulltext"
            Retrieve the interactions in the article from the source (fulltext or abstract).
        '''
        self.extract_sentences(source=source)
        for sentence in self.sentences:
            sentence.annotate()
            sentence.get_candidates(only_dict)
            for candidate in sentence.candidates:
                candidate.predict()

    def as_html(self):
        '''
        Writes tokenized sentences as HTML
        '''
        pass

    def extract_sentences(self, mode="split", source="fulltext"):
        '''
        Finds sentence boundaries and saves them as sentence objects in the
        attribute "sentences" as a list of Sentence objects.

        Parameters
        ----------
        mode : str, optional, default = "split"
            Split the sentences ("split") or use the whole "source" as a single sentence ("no-split").
            Useful for developing and debugging.

        source : str, optional, default = "fulltext"
            Use the "fulltext" or the "abstract" to extract sentences.
        '''
        text = ""
        if source == "fulltext":
            text = str(self.fulltext)
        else:
            text = str(self.abstract)

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
            sentences = [ s.strip() for s in sentences ]
            try:
                html = HTMLParser()
            except:
                import html
            for sentence in sentences:
                sentence = str(html.unescape(sentence))
                if not sentence.strip() or not isinstance(sentence, str):
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

    def __str__(self):
        return "Article with PubMED id:%s" % (self.pmid)

# ----------------------------------------------
class Protein(object):
    '''
    Class for protein in sentences.

    Attributes
    ----------
    symbol : str, no default
        Symbol of the protein or the gene.

    positions : list, no default
        List of the position of the protein in the tokenized sentence (1-Indexed).

    sentence : Sentence, no default
        Sentence object to which the protein belongs.

    synonym : str, no default
        Synonymous symbol of the protein/gene.

    count : int, no default
        Length of position list.

    '''
    GENEDICT = dict()
    GENEDICTFILE = pkg_resources.resource_filename('ppaxe', 'data/HGNC_gene_dictionary.txt')
    try:
        with open(GENEDICTFILE, 'r') as f:
            for line in f:
                line = line.strip()
                cols = line.split("\t")
                GENEDICT[cols[0]] = cols[0]
                for alias in cols[1:]:
                    GENEDICT[alias.upper()] = cols[0]
    except Exception:
        raise GeneDictError("Can't read %s\n" % GENEDICTFILE)


    def __init__(self, symbol, positions, sentence):
        '''
        Parameters
        ----------
        symbol : str, required, no default
            Symbol of the protein or the gene.

        positions : list, required, no default
            List of the position of the protein in the tokenized sentence (1-Indexed).

        sentence : Sentence, required, no default
            Sentence object to which the protein belongs.
        '''

        #self.disambiguate()
        self.symbol = symbol
        self.positions = positions
        self.sentence = sentence
        self.synonym = list()
        self.count = len(positions)

    def is_in_dict(self):
        '''
        Checks if a given protein symbol is in dictionary
        '''
        disambiguated = self.symbol.upper()
        disambiguated = disambiguated.replace("'", "")
        disambiguated = disambiguated.replace('"', '')
        if disambiguated in Protein.GENEDICT:
            return True
        else:
            return False

    def disambiguate(self):
        '''
        Method for disambiguating the gene (convert it to the approved symbol if possible).
        '''
        disambiguated = self.symbol.upper()
        disambiguated = disambiguated.replace("'", "")
        disambiguated = disambiguated.replace('"', '')
        if disambiguated in Protein.GENEDICT:
            return Protein.GENEDICT[disambiguated]
        else:
            return disambiguated

    def __str__(self):
        return "%s found in positions %s" % (self.symbol, ":".join([ str(idx) for idx in self.positions ]))


# ----------------------------------------------
class Sentence(object):
    '''
    Class for sentences.

    Attributes
    ----------
    originaltext : str, no default
        Original text string of the sentence.

    tokens : list, no default
        List of tokens retrieved from StanfordCoreNLP. Each element is a dictionary
        with keys:
            "index" : Position of token (1-Indexed).
            "word"  : Word of the token.
            "lemma" : Lemma of the token.
            "ner"   : Protein ("P") or Other ("O").
            "pos"   : Part-of-Speech tag.

    candidates : list, no default
        List of Candidate objects in sentence.

    proteins : list, no default
        List of Protein objects found in sentence.

    '''
    def __init__(self, originaltext):
        '''
        Parameters
        ----------
        originaltext : str, required, no default
            Original text string of the sentence.
        '''
        self.originaltext = originaltext
        self.tokens       = list()
        self.tree         = list()
        self.candidates   = list()
        self.proteins     = list()

    def annotate(self):
        '''
        Annotates the genes/proteins in the sentence using StanfordCoreNLP
        trained NER tagger. Will add a list of tokens to the attribute "tokens".
        '''
        if not self.originaltext.strip():
            self.tokens = ""
        annotated = json.loads(NLP.annotate(self.originaltext))
        if annotated['sentences']:
            self.tokens = annotated['sentences'][0]['tokens']

    def get_candidates(self, only_dict=False):
        '''
        Gets interaction candidates candidates for sentence (attribute: candidates)
        and all the proteins (attribute: proteins).
        '''
        if not self.tokens:
            self.annotate()
        if self.candidates:
            return
            
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
            if only_dict is True:
                if not protein.is_in_dict():
                    continue
            self.proteins.append(protein)
            prots_in_sentence.append(protein)

        # Create candidates for sentence
        for prot in itertools.combinations(prots_in_sentence, r=2):
            self.candidates.append(InteractionCandidate(prot1=prot[0], prot2=prot[1]))

    def to_html(self):
        '''
        Sentence to HTML string tagging the proteins and the verbs using <span> tags.
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
    FEATURES: in feature_names.py.

    Attributes
    ----------
    prot1 : Protein, no default
        Protein object of the first protein involved in the possible interaction.

    prot2 : Protein, no default
        Protein object of the second protein involved in the possible interaction.

    between_idxes : tuple, no default
        Indexes of end of Protein_1 and start of Protein_2 (1-Indexed).

    label : bool, no default
        Label of Candidate when prediction is performed. True for interacting proteins
        and False for non-interacting proteins. True if votes >= 0.55.

    votes : float, no default
        Percentage of votes of the Random Forest Classifier.

    feat_cols : list, no default
        Feature column indexes of the non-zero features computed for Candidate.

    feat_current_col : int, no default
        Store the current feature column index that has been computed.

    feat_vals : list, no default
        Values of the non-zero features.

    features_sparse : sparse.coo_matrix, no default
        Sparse Coo matrix with features for Candidate.

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
        try:
            predictor = pickle.load(f)
        except:
            try:
                predictor = pickle.load(f, encoding='latin1')
            except Exception as err:
                raise(CannotLoadClassifier("RF_scikit.pkl can't be loaded - %s" % err))



    def __init__(self, prot1, prot2):
        '''
        Parameters
        ----------
        prot1 : Protein, required, no default
            Protein object of the first protein involved in the possible interaction.

        prot2 : Protein, required, no default
            Protein object of the second protein involved in the possible interaction.
        '''
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
        is a real interaction. Fills attribute features_sparse, which is a Scipy sparse matrix.
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
        Total tokens in sentence.
        '''
        if len(self.prot1.sentence.tokens) > 0:
            self.feat_cols.append(self.feat_current_col)
            self.feat_vals.append(len(self.prot1.sentence.tokens))
        self.feat_current_col += 1

    def __prot_count(self, mode="all"):
        '''
        Counts the number of times the proteins of the candidate (A and B)
        appear in the sentence (either in the whole sentence [mode="all"] or between
        A and B [mode="between"] ).

        Parameters
        ----------
        mode : str, optional, default = "all"
            Count number of times proteins appears in whole sentence (mode="all") or only
            between candidate proteins (mode="between").
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
        Returns a string with the token POS annotations for coordinates 'from'-'to'.

        Parameters
        ----------
        mode : str, optional, default = "all"
            Protein object of the first protein involved in the possible interaction.

        prot2 : Protein, required, no default.
            Protein object of the second protein involved in the possible interaction.

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

    def __pos_features(self, mode):
        '''
        Gets counts of POS tags and adds those features

        Parameters
        ----------
        mode : str, required, no default
            Counts POS tags in whole sentence (mode="all") or between candidate proteins (mode="between").
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
        Computes the two verb distances (closest and farthest).

        Parameters
        ----------
        pidx : int, required, no default
            Index of the protein in tokenized sentence (1-Indexed).

        vidxes : list, required, no default
            List of indexes (ints) of the verbs in tokenized sentence (1-Indexed).
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
        and for only the verbs between protein 1 and protein 2 (flag=between).

        Parameters
        ----------
        flag : str, required, no default
            Compute verb features for whole sentence (flag="all") or only between
            candidate proteins (flag="between").
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
        Counts appearance of keywords in sentence.

        Parameters
        ----------
        mode : str, optional, default = "all"
            Count keywords in whole sentence (mode="all") or between candidate proteins (mode="between").
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

        for word, value in sorted(keywords.items()):
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
        Computes the votes (prediction) of the candidate by using the Random Forest
         classifier trained with scikitlearn.
        '''
        if self.features_sparse is None:
            self.compute_features()
        pred = InteractionCandidate.predictor.predict_proba(self.features_sparse)[:,1]

        self.votes = pred[0]
        self.__normalize_pred()
        if self.votes >= 0:
            self.label = True
        else:
            self.label = False


    def __normalize_pred(self):
        '''
        Changes scale of votes in self.votes from 0.55-1 to 0-1.
        '''
        before = self.votes
        self.votes = (self.votes-0.55)/(1-0.55)
        self.votes = round(self.votes, 3)

        


    def to_html(self):
        '''
        Transforms candidate to html with only involved proteins tagged and only
        verbs between proteins tagged.
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

class GeneDictError(Exception):
    '''
    Raised when dictionary can't be read
    '''
    pass

class CannotLoadClassifier(Exception):
    '''
    Raised when classifier can't be loaded
    '''
    pass