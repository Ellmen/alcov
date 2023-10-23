from setuptools import setup

with open('README.md', 'r', encoding='utf-8') as fh:
    long_description = fh.read()

setup(
    name='alcov',
    version='1.1.19',
    description='Identify frequencies of concerning mutations from aligned reads',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Isaac Ellmen',
    author_email='isaac.ellmen@uwaterloo.ca',
    maintainer='Jenn Knapp',
    maintainer_email='jenn.knapp@uwaterloo.ca',
    packages=['alcov'],
    exclude_package_data={
	'alcov': ['data'],
	},
    url='https://github.com/Ellmen/alcov',
    install_requires=[
        'fire',
        'numpy',
        'pandas',
        'scikit-learn>=0.24',
        'matplotlib',
        'seaborn',
        'pysam',
    ],
    entry_points={
        'console_scripts': ['alcov=alcov.command_line:main'],
    }
)
