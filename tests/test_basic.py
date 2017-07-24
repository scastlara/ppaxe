# -*- coding: utf-8 -*-
'''
Tests for the main classes of ppaxe
'''

from ppaxe import core


def test_single_article_query():
    '''
    Test if single article query for fulltext in PMC works
    '''
    query = core.PMQuery(ids=["4304705"], database="PMC")
    query.get_articles()
    assert(query.articles[0].pmid == "25615823")

def test_multiple_article_query():
    '''
    Test if multiple queries work!
    '''
    query = core.PMQuery(ids=["4304705","5055395"], database="PMC")
    query.get_articles()
    pmid_concatenation = query.articles[0].pmid + query.articles[1].pmid
    assert(pmid_concatenation == "2561582327612382")
