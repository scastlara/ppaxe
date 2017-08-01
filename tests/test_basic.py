# -*- coding: utf-8 -*-
'''
Tests for the main classes of ppaxe
'''
from ppaxe import core
from pycorenlp import StanfordCoreNLP
import json

def test_sentence_separator():
    '''
    Tests if sentence separator works...
    '''
    article = core.Article(pmid="1234", fulltext="""
        To identify roles of Hh signaling in the planarian CNS maintenance, we examined gene expression changes using RNA sequencing of cephalic ganglia following RNAi of hh, ptc, or a control gene (C. elegans unc-22) not present in the planarian genome.
        We developed a dissection technique that allowed cephalic ganglia tissue to be collected from large (>2 cm) S2F1L3F2 sexual strain S. mediterranea animals following a brief acid-based fixation (Figure 1C).
        To test for enrichment using this dissection technique, amputated head fragments collected from CIW4 asexual strain S. mediterranea animals after six control dsRNA feedings were used as a reference library (Figure 1D).
        Head fragments contain cephalic ganglia as well as most major planarian tissue types (Hyman, 1951).
        The magic number is 12.45 for the species S. mediterranea.
        Figure 2.a and 3.B Is the most important.
        S. mediterranea and C. elegans.
        But not (S.mediterranea)
    """)
    article.extract_sentences()
    assert(len(article.sentences) == 8)

def test_stanford_corenlp_server():
    '''
    Tests connection to stanford corenlp server
    '''
    try:
        nlp = StanfordCoreNLP('http://localhost:9000')
        nlp.annotate("HOLA")
        assert(1 == 1)
    except:
        assert("Connection error" == "StanfordCoreNLP")

def test_annotation():
    '''
    Tests if Stanford coreNLP NER works
    '''
    text = "MAPK seems to interact with chloroacetate esterase"
    nlp = StanfordCoreNLP('http://localhost:9000')
    ner_list = list()
    for token in json.loads(nlp.annotate(text))['sentences'][0]['tokens']:
        ner_list.append(token['ner'])
    assert("".join(ner_list) == "POOOOPP")

def test_article_annotation():
    '''
    Tests if article annotation works
    '''
    article_text = """
        MAPK seems to interact with chloroacetate esterase.
        However, MAPK is a better target for peroxydase.
        The thing is, Schmidtea mediterranea is a good model organism because reasons.
        However, cryoglobulin is better.
    """
    prot_list = list()
    article = core.Article(pmid="1234", fulltext=article_text)
    article.annotate_sentences()
    for sentence in article.annotated:
        for token in sentence:
            if token['ner'] == "P":
                prot_list.append(token['word'])
    assert(",".join(prot_list) == "MAPK,chloroacetate,esterase,MAPK,peroxydase,cryoglobulin")

def test_get_proteins():
    '''
    Tests the retrieval of candidates
    '''
    text = "However, MAPK is a better target for chloroacetate esterase which is an essential protein for cryoglobulin."
    article = core.Article(pmid="1234", fulltext=text)
    article.annotate_sentences()
    article.get_candidates()
    assert(str(article.candidates[0]) == "[MAPK] may interact with [chloroacetate esterase]")

def test_token_distance():
    '''
    Tests token distance feature
    '''
    text = "The protein MAPK interacts directly with cryoglobulin which is very interesting."
    article = core.Article(pmid="1234", fulltext=text)
    article.annotate_sentences()
    article.get_candidates()
    for candidate in article.candidates:
        candidate.compute_features()
        assert(candidate.features[0] == 3)

def test_total_tokens():
    '''
    Tests total tokens
    '''
    text = "The protein MAPK interacts directly with cryoglobulin which is very interesting."
    article = core.Article(pmid="1234", fulltext=text)
    article.annotate_sentences()
    article.get_candidates()
    for candidate in article.candidates:
        candidate.compute_features()
        assert(candidate.features[1] == 12)

def test_verb_features():
    '''
    Tests total number of verbs between proteins
    '''
    text = "The protein MAPK is interacting and activating directly with cryoglobulin which is very interesting."
    article = core.Article(pmid="1234", fulltext=text)
    article.annotate_sentences()
    article.get_candidates()
    for candidate in article.candidates:
        candidate.compute_features()
        # Check if it has detected one VBG (interacting) verb and another VBZ (is) verb
        # and verb scores (for now 3 and 5)
        assert(candidate.features[2:14] == [0, 0, 2, 0, 0, 1, 4, 7, 1, 4, 3, 6])

def test_candidates_multiple_sentences():
    '''
    Tests if multiple sentences in one string work
    '''
    text = "In patients with a complete response to PROT4  , the levels of PROT2  were higher at 24 weeks following PROT4  treatment than that of pre - treatment ( P = 0.04 ) , and the levels of PROT3  decreased markedly at 12 and 24 weeks ( P = 0.02 , 0.03 , respectively ) . mRNA expression positively correlated with the level of PROT55 / Th2 type cytokines in the PROT99 ."
    article = core.Article(pmid="1234", fulltext=text)
    article.annotate_sentences()
    article.get_candidates()
    assert(article.candidates[-1].prot2.symbol == "PROT99")


def test_get_pos_annotation():
    '''
    Tests the method get_token_pos()
    '''
    text = "The protein MAPK14 seems to interact with MAPK12."
    article = core.Article(pmid="1234", fulltext=text)
    article.annotate_sentences()
    article.get_candidates()
    assert(article.candidates[0]._InteractionCandidate__get_token_pos(mode="between") == "NN,VBZ,TO,VB,IN,NN")

def get_pos_count():
    '''
    Tests pos count feature
    '''
    text = "The protein MAPK14 seems to interact with MAPK12."
    article = core.Article(pmid="1234", fulltext=text)
    article.annotate_sentences()
    article.get_candidates()
    for candidate in article.candidates:
        candidate.compute_features()
        assert(candidate.features[85] == 1)
