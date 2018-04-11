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
    pmid_concatenation = set([query.articles[0].pmid, query.articles[1].pmid])
    assert(pmid_concatenation == set(["25615823","27612382"]))

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

def test_year_pmc():
    '''
    Tests retrieval of year from PMC XML response
    '''
    query = core.PMQuery(ids=["25615823"], database="PMC")
    query.get_articles()
    for article in query:
        assert(article.year == "2015")

def test_year_pubmed():
    '''
    Tests retrieval of year from PMC XML response
    '''
    query = core.PMQuery(ids=["25615823"], database="PUBMED")
    query.get_articles()
    for article in query:
        assert(article.year == "2015")

def test_article_journal_pmc():
    '''
    Tests the retrieval of the article journal from the PMC XML response
    '''
    query = core.PMQuery(ids=["25615823"], database="PMC")
    query.get_articles()
    for article in query:
        assert(article.journal == "PLoS One")


def test_article_journal_pubmed():
    '''
    Tests the retrieval of the article journal from the PubMed XML response
    '''
    query = core.PMQuery(ids=["28869924"], database="PUBMED")
    query.get_articles()
    for article in query:
        assert(article.journal == "Aquatic toxicology (Amsterdam, Netherlands)")

def test_too_many_pmids():
    '''
    Tests if ppaxe can handle many PMID in pmid_2_pmc
    '''
    identifiers = ["26267445"] * 300
    pmcids = core.pmid_2_pmc(identifiers)
    assert(len(pmcids) == 1)

def test_too_many_pmids_2():
    '''
    Tests if ppaxe can handle a lot of PMIDs (more than GET can handle)
    '''
    identifiers = ["26267445"] * 2000
    query = core.PMQuery(ids=identifiers, database="PUBMED")
    query.get_articles()
    assert(query.found == set(["26267445"]))
