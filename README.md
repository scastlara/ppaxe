
<img width="350" src="https://raw.githubusercontent.com/scastlara/ppaxe/master/ppaxe/logo.png"/>

-----

Tool to retrieve **protein-protein interactions** and calculate protein/gene symbol ocurrence in the scientific literature (PubMed & PubMedCentral).


## Prerequisites


```
xml.dom
numpy
pycorenlp
cPickle
scipy
```

## Installing

#### Install ppaxe
You can install this package using _pip_. However, before doing so, you have to download the [Random Forest predictor](https://www.dropbox.com/s/t6qcl19g536c0zu/RF_scikit.pkl?dl=0) and place it in `ppaxe/data`.

```
# Clone the repository
git clone https://github.com/scastlara/ppaxe.git

# Download pickle with RF
wget https://www.dropbox.com/s/t6qcl19g536c0zu/RF_scikit.pkl?dl=0 -O ppaxe/ppaxe/data/RF_scikit.pkl

# Install
pip install ppaxe
```

#### Download StanfordCoreNLP
In order to use the package you will need a [StanfordCoreNLP](https://stanfordnlp.github.io/CoreNLP) server setup with
 the [Protein/gene Tagger](https://www.dropbox.com/s/ec3a4ey7s0k6qgy/FINAL-ner-model.AImed%2BMedTag%2BBioInfer.ser.gz?dl=0).

 ```
 # Download StanfordCoreNLP
 wget http://nlp.stanford.edu/software/stanford-corenlp-full-2017-06-09.zip
 unzip stanford-corenlp-full-2017-06-09.zip

 # Download the Protein tagger
 wget https://www.dropbox.com/s/ec3a4ey7s0k6qgy/FINAL-ner-model.AImed%2BMedTag%2BBioInfer.ser.gz?dl=0 -O FINAL-ner-model.AImed+MedTag+BioInfer.ser.gz

 # Download English tagger models
 wget http://nlp.stanford.edu/software/stanford-english-corenlp-2017-06-09-models.jar -O stanford-corenlp-full-2017-06-09/stanford-english-corenlp-2017-06-09-models.jar

 # Change the location of the tagger in ppaxe/data/server.properties if necessary
 # ...

 # Start the StanfordCoreNLP server
 cd stanford-corenlp-full-2017-06-09/
java -mx1000m -cp ./stanford-corenlp-3.8.0.jar:stanford-english-corenlp-2017-06-09-models.jar edu.stanford.nlp.pipeline.StanfordCoreNLPServer -port 9000 -serverProperties ~/ppaxe/ppaxe/data/server.properties
 ```

Once the server is up and running and ppaxe has been installed, you are good to go.

By default, ppaxe will assume the server is available at localhost:9000. If you want to change the address, set up the server with the appropiate port and change the address in ppaxe by assigning the new address to the variable ppaxe.ppcore.NLP:

* **Start the server**

```
# Change the location of the ner tagger in server.properties manually
java -mx10000m -cp ./stanford-corenlp-3.8.0.jar:stanford-english-corenlp-2017-06-09-models.jar edu.stanford.nlp.pipeline.StanfordCoreNLPServer -port your_port -serverProperties ppaxe/data/server.properties
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
ppaxe -p pmids.txt -d PMC -v -o output.tbl -r report
```

### Report

The report output (`option -r`) will contain a simple summary of the analysis, the interactions retrieved (including the sentences from which they were retrieved), a table with the protein/gene counts and a graph visualization made using [cytoscape.js](http://js.cytoscape.org/).

<img src="https://raw.githubusercontent.com/scastlara/ppaxe/master/ppaxe/data/report1-example.png"/>
<img src="https://raw.githubusercontent.com/scastlara/ppaxe/master/ppaxe/data/report2-example.png"/>

## Documentation

Refer to the [wiki](https://github.com/scastlara/ppaxe/wiki/Documentation) of the package.

## Running the tests

To run the tests:

```
python -m pytest -v tests
```

## Authors

* **Sergio Castillo-Lara** - at the [Computational Genomics Lab](https://compgen.bio.ub.edu)


## License

This project is licensed under the GNU GPL3 license - see the [LICENSE](LICENSE) file for details
