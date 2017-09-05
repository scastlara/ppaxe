
<img width="350" src="https://raw.githubusercontent.com/scastlara/ppaxe/master/ppaxe/logo.png"/>

-----

Tool to retrieve **protein-protein interactions** and calculate protein/gene symbol ocurrence in the scientific literature (PubMed & PubMedCentral).


## Prerequisites


```
xml.dom
numpy
pycorenlp
markdown
cPickle
scipy
```

## Installing

You can install this package using _pip_

```
git clone https://github.com/scastlara/ppaxe.git
pip install ppaxe
```

In order to use the package you will need a [StanfordCoreNLP](https://stanfordnlp.github.io/CoreNLP) server setup with
 the [Protein/gene Tagger](https://compgen.bio.ub.edu). By default, it will assume it is available at localhost:9000. If you want to change the address you can do it as follows:

* **Start the server**

```
# Change the location of the ner tagger in server.properties manually
java -mx10000m -cp ./stanford-corenlp-3.8.0.jar:stanford-english-corenlp-2017-06-09-models.jar edu.stanford.nlp.pipeline.StanfordCoreNLPServer -port 9000 -serverProperties ppaxe/data/server.properties
```

* **Use the ppaxe package**


 ```py

from ppaxe import core as ppcore
from pycorenlp import StanfordCoreNLP

ppcore.NLP = StanfordCoreNLP(your_new_adress)

# Do whatever you want
 ```

## Usage

### ppaxe core classes

```py
from ppaxe import core as ppcore

# Perform query to PubMedCentral
pmids = ["28615517","28839427","28831451","28824332","28819371","28819357"]
query = ppcore.PMQuery(ids=pmids, database="PMC")
query.get_articles()

# Retrieve interactions from text
for article in query:
    article.predict_interactions()

# Iterate through predictions
for article in query:
    for sentence in article.sentences:
        for candidate in sentence.candidates:
            if candidate.label is True:
                # We have an interaction
                print("%s interacts with %s in article %s" % (candidate.prot1.symbol, candidate.prot2.symbol, article.pmid ))
                print(candidate.to_html())

# Print html report
# Will create 'report_file.html'
summary = ppcore.ReportSummary(query)
summary.make_report("report_file")
```

### ppaxe script

```sh
# Will read PubMed ids in pmids.txt, predict the interactions
# in their fulltext from PubMedCentral, and print a tabular output
# and an html report
bin/ppaxe -p pmids.txt -d PMC -v -o output.tbl -r report
```

## Documentation

Refer to the [wiki](https://github.com/scastlara/ppaxe/wiki) of the package.

## Running the tests

To run the tests:

```
python -m pytest -v tests
```

## Authors

* **Sergio Castillo-Lara** - at the [Computational Genomics Lab](https://compgen.bio.ub.edu)


## License

This project is licensed under the GNU GPL3 license - see the [LICENSE](LICENSE) file for details
