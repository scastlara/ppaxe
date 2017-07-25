# -*- coding: utf-8 -*-
'''
Test connections
'''
from ppaxe import core

def test_single_article_query():
    '''
    Test if single article query for fulltext in PMC works
    '''
    query = core.PMQuery(ids=["25615823"], database="PMC")
    query.get_articles()
    assert(query.articles[0].pmid == "25615823")

def test_multiple_article_query():
    '''
    Test if multiple queries work!
    '''
    query = core.PMQuery(ids=["25615823","27612382"], database="PMC")
    query.get_articles()
    pmid_concatenation = query.articles[0].pmid + query.articles[1].pmid
    assert(pmid_concatenation == "2561582327612382")

def test_query_text():
    '''
    Tests if fulltext is correctly retrieved
    '''
    query = core.PMQuery(ids=["25615823"], database="PMC")
    query.get_articles()
    known_text = "Liver cancer is the sixth most frequent"
    assert(query[0].fulltext.split("\n")[0].split(".")[0][0:39] == known_text)

def test_query_found():
    '''
    Tests if PMQuery separates input ids as found and not found
    '''
    query = core.PMQuery(ids=["99999999","27612382"], database="PMC")
    query.get_articles()
    assert(query.found == set(["27612382"]))

def test_query_notfound():
    '''
    Tests if PMQuery separates input ids as found and not found
    '''
    query = core.PMQuery(ids=["99999999","27612382"], database="PMC")
    query.get_articles()
    assert(query.notfound == set(["99999999"]))
