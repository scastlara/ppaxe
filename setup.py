import setuptools

requires = [
    'pycorenlp',
    'scipy',
    'sklearn',
    'requests',
    'xml',
    'networkx',
    'uuid',
    'matplotlib'
    ]

setuptools.setup(name='ppaxe',
      version='1.0',
      description='Protein-Protein interactions extractor from PubMed articles',
      url='http://github.com/scastlara/ppaxe',
      author='S. Castillo-Lara',
      author_email='s.cast.lara@gmail.com',
      license='GPL-3.0',
      scripts=['bin/ppaxe'],
      include_package_data=True,
      packages=setuptools.find_packages(),
      install_requires=requires,
      package_data = { 'ppaxe' : ['data/RF_scikit.pkl', 'data/HGNC_gene_dictionary.txt', 'data/cytoscape_template.js', 'data/style.css']},
      zip_safe=False)

