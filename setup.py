from setuptools import setup

with open('README.md', 'r', encoding='utf-8') as fh:
    long_description = fh.read()

setup(
    name='alcov',
    version='0.1.1',
    description='Identify frequencies of concerning mutations from aligned reads',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Isaac Ellmen',
    author_email='isaac.ellmen@uwaterloo.ca',
    packages=['alcov'],
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
