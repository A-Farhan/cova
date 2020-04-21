import setuptools

def readme():
    with open('README.rst') as f:
        return f.read()

setuptools.setup(name='cova',
        version='0.0.1',
        license='MIT License',
        author='Farhan Ali',
        author_email='afarhan@ncbs.res.in',
        description='coronavirus variant analysis (simplified)',
        keywords='variants diversity selection phylogeny coronavirus covid19',
        long_description=readme(),
        url='https://bitbucket.org/Farhan_Bugbears/cova',
        packages=setuptools.find_packages(),
        classifiers=[
            "Programming Language :: Python :: 3",
            "License :: OSI Approved :: MIT License",
            "Operating System :: Linux",
        ],
        python_requires='>=3.6',
        install_requires=[ 
        	'biopython>=1.75',
        	'click',
        	'matplotlib>=3.1',
        	'numpy>=1.17',
        	'scipy>=1.3',
        	'seaborn',
        	'ete3>=3.1.1',
            'PyQt5',
        	'joblib>=0.14.0'
        ],
        test_suite='nose.collector',
    	tests_require=['nose'],
        package_dir={'cova':'cova'},
        package_data={'cova':['data/','testdata/', 'softwares/']},
        include_package_data=True,
        zip_safe=False,
        entry_points = {
        'console_scripts': ['CoVa=cova.command_line:cli'],
    }
    )
