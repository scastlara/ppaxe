# -*- coding: utf-8 -*-
'''
Tests for the main classes of ppaxe
'''
from ppaxe import core


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
