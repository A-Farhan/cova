import setuptools

def readme():
    with open('README.rst') as f:
        return f.read()

setuptools.setup(name='cova',
        version='0.3',
        license='MIT License',
        author='Farhan Ali',
        author_email='afarhan@ncbs.res.in',
        description='coronavirus variant analysis (simplified)',
        keywords='variants diversity selection phylogeny coronavirus covid19',
        long_description=readme(),
        url='https://github.com/A-Farhan/cova',
        packages=setuptools.find_packages(),
        classifiers=[
            "Programming Language :: Python :: 3",
            "License :: OSI Approved :: MIT License",
            "Operating System :: Linux",
        ],
        python_requires='>=3.6',
        install_requires=[ 
        	'biopython>=1.75',
            'pandas>=1.0.5'
        	'click',
        	'matplotlib>=3.1',
        	'numpy>=1.17',
        	'scipy>=1.3',
            'seaborn', 
        ],
        extras_require={
            'treeplot':[ 'ete3>=3.1.1', 'PyQt5'],
        },
        test_suite='nose.collector',
    	tests_require=['nose'],
        package_dir={'cova':'cova'},
        package_data={'cova':['data/','testdata/']},
        include_package_data=True,
        zip_safe=False,
        entry_points = {
        'console_scripts': [
            'CoVa=cova.command_line:cli', 
            'preprocess=cova.preprocessing:main_fun',
            'replace_gisaid_header=cova.replace_fasta_header_gisaid_accession:main_fun',
            'extract_metadata=cova.extract_metadata:main_fun',
            'plottree=cova.plottree:main_fun [treeplot]'
        ],
    }
)
