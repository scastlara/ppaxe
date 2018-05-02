import setuptools

setuptools.setup(name='ppaxe',
      version='0.2',
      description='PPI extractor from PubMed articles',
      url='http://github.com/scastlara/ppaxe',
      author='S. Castillo-Lara',
      author_email='s.cast.lara@gmail.com',
      license='GPL-3.0',
      scripts=['bin/ppaxe'],
      include_package_data=True,
      packages=setuptools.find_packages(),
      package_data = { 'ppaxe' : ['data/RF_scikit.pkl', 'data/HGNC_gene_dictionary.txt', 'data/cytoscape_template.js', 'data/style.css']},
      zip_safe=False)

