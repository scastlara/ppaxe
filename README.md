
<img width="350" src="https://raw.githubusercontent.com/scastlara/ppaxe/master/ppaxe/logo.png"/>

-----

# PP-axe

Tool to retrieve protein-protein interactions and calculate protein/gene symbol ocurrence in the scientific literature (PubMed & PubMedCentral).


### Prerequisites


```
xml.dom
numpy
pycorenlp
```

### Installing

You can install this package using _pip_

```
git clone https://github.com/scastlara/ppaxe.git
pip install ppaxe
```

In order to use the package you will need a (StanfordCoreNLP)[https://stanfordnlp.github.io/CoreNLP/] server setup with
 the [Protein/gene Tagger](https://compgen.bio.ub.edu). By default, it will assume it is available at localhost:9000. If you want to change the address you can do it as follows:

* **Start the server**

```
# Change the location of the ner tagger in server.properties manually
java -mx10000m -cp ./stanford-corenlp-3.8.0.jar:stanford-english-corenlp-2017-06-09-models.jar edu.stanford.nlp.pipeline.StanfordCoreNLPServer -port 9000 -serverProperties ppaxe/data/server.properties
```

* **Use the ppaxe package**
 ```py
from ppaxe import core
from pycorenlp import StanfordCoreNLP

core.NLP = StanfordCoreNLP(your_new_adress)

# Do whatever you want
 ```

## Running the tests

To run the tests:

```
python -m pytest ppaxe/tests
```

## Authors

* **Sergio Castillo-Lara** - at the [Computational Genomics Lab](https://compgen.bio.ub.edu)


## License

This project is licensed under the GNU GPL3 license - see the [LICENSE](LICENSE) file for details
